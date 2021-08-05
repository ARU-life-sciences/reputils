pub mod dot {

    use bio::io::fasta;
    use clap::value_t;
    use plotters::prelude::*;

    use crate::utils::windows::SeqWindows;

    // the code here is based on the laconic R version here:
    // https://github.com/cran/seqinr/blob/master/R/dotPlot.R
    // It's not an optimal algorithm by any stretch

    pub fn dot(matches: &clap::ArgMatches) -> Result<(), Box<dyn std::error::Error + 'static>> {
        // load in the fasta
        let fasta = matches.value_of("fasta").unwrap();
        let wsize = value_t!(matches.value_of("wsize"), usize).unwrap_or_else(|e| e.exit());
        let wstep = value_t!(matches.value_of("wstep"), usize).unwrap_or_else(|e| e.exit());
        let nmatch = value_t!(matches.value_of("nmatches"), usize).unwrap_or_else(|e| e.exit());
        let dir = value_t!(matches.value_of("dir"), String).unwrap_or_else(|e| e.exit());

        let reader = fasta::Reader::from_file(fasta).expect("[-]\tPath invalid.");

        for record in reader.records() {
            let record = record.expect("[-]\tError during fasta record parsing.");
            let seq_counter = SeqWindows::new(record.seq(), wsize, wstep);
            let no_its = seq_counter.count();

            let seq_windows = SeqWindows::new(record.seq(), wsize, wstep);

            let seq_windows_vec: Vec<&[u8]> = seq_windows.map(|s| s).collect();

            // Base 1d array
            let mut matrix_raw = vec![false; no_its * no_its];
            // Vector of 'width' elements slices
            let mut matrix_base: Vec<_> = matrix_raw.as_mut_slice().chunks_mut(no_its).collect();
            // Final 2d array `&mut [&mut [_]]`
            let matrix: &mut [&mut [bool]] = matrix_base.as_mut_slice();

            // compare windows to themselves
            for row in 0..no_its {
                for column in 0..no_its {
                    // base or kmer 1
                    let k1 = seq_windows_vec[row];
                    // base or kmer 2
                    let k2 = seq_windows_vec[column];
                    matrix[row][column] = match_case(k1, k2, nmatch)
                }
            }
            matplot(matrix, &dir, record.id())?;
        }
        Ok(())
    }

    fn matplot(
        matrix: &mut [&mut [bool]],
        dir: &str,
        name: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let path = format!("{}/{}.png", dir, name);
        let dim = 1024u32;
        let root = BitMapBackend::new(&path, (dim, dim)).into_drawing_area();

        root.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(&root)
            .margin(5)
            .top_x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(0i32..(matrix.len() as i32), (matrix.len() as i32)..0i32)?;

        chart
            .configure_mesh()
            .x_labels(15)
            .y_labels(15)
            .x_label_offset(35)
            .y_label_offset(25)
            .disable_x_mesh()
            .disable_y_mesh()
            .label_style(("sans-serif", 20))
            .draw()?;

        chart.draw_series(
            matrix
                .iter()
                .zip(0..)
                .map(|(l, y)| l.iter().zip(0..).map(move |(v, x)| (x as i32, y as i32, v)))
                .flatten()
                .map(|(x, y, v)| {
                    Rectangle::new(
                        [(x, y), (x + 1, y + 1)],
                        RGBColor(
                            match v {
                                true => 0u8,
                                false => 255u8,
                            },
                            match v {
                                true => 0u8,
                                false => 255u8,
                            },
                            match v {
                                true => 0u8,
                                false => 255u8,
                            },
                        )
                        .filled(),
                    )
                }),
        )?;

        Ok(())
    }

    fn match_case<'a>(kmer1: &'a [u8], kmer2: &'a [u8], nmatch: usize) -> bool {
        // if kmer1 or kmer2 contains a N or - skip?
        if kmer1.contains(&45u8)
            || kmer2.contains(&45u8)
            // little n
            || kmer1.contains(&110u8)
            || kmer2.contains(&110u8)
            // big N
            || kmer1.contains(&78u8)
            || kmer2.contains(&78u8)
        {
            return false;
        }

        // iterate over both kmers at the same time
        let mut counts = Vec::new();
        for (b1, b2) in kmer1.iter().zip(kmer2.iter()) {
            if b1 == b2 {
                counts.push(0)
            } else {
                counts.push(1)
            }
        }
        let counts_sum: usize = counts.iter().sum();
        // if sum less than or equal to number of matches,
        // we want to return true
        counts_sum <= nmatch
    }
}
