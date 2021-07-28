pub mod div {
    use bio::io::fasta;
    use clap::value_t;
    use plotters::prelude::*;

    // heavily poached from https://github.com/noahaus/sliding-window-scripts/blob/13d872b379a3501ca9b506f8bcccf89a9cd81c8d/tjd/src/main.rs

    pub fn diversity_windows(matches: &clap::ArgMatches) {
        let fasta = matches.value_of("fasta").unwrap();
        let window_size = value_t!(matches.value_of("window"), usize).unwrap_or_else(|e| e.exit());
        let window_step = value_t!(matches.value_of("step"), usize).unwrap_or_else(|e| e.exit());
        let plot = matches.is_present("plot");

        let mut reader = fasta::Reader::from_file(fasta)
            .expect("[-]\tPath invalid.")
            .records();

        let mut alignment = Vec::<String>::new();

        while let Some(Ok(record)) = reader.next() {
            alignment.push(std::str::from_utf8(record.seq()).unwrap().to_string());
        }

        let mut window = Vec::<String>::new(); // vector to store window
        let mut start: usize = 0; // index to start the window
        let mut end: usize = start + window_size; // index to end the window

        let mut data = Vec::new();

        while end <= alignment[0].len() {
            for row in &alignment {
                window.push(row[start..end].to_string());
            }

            let pi = calculate_pi(&window);
            // window contains a vec of chunks.
            println!("{}\t{}\t{:.3}", start, end, pi);
            data.push((start, end, pi));

            start += window_step;
            end += window_step;
            window.clear();
        }

        if plot {
            // do plot
            plot_pi_windows(data).unwrap();
        }
    }

    fn calculate_pi(alignment: &Vec<String>) -> f32 {
        // number of sequences being analyzed
        let align_length = alignment.len();
        // vector that holds all the values of the pairwise distances
        let mut distances = Vec::new();

        // calculate hamming distance between two sequences at a time
        for i in 0..align_length {
            for j in 0..align_length {
                if i < j {
                    distances.push(hamming_distance(&alignment[i], &alignment[j]));
                } else {
                    continue;
                }
            }
        }

        // sum the pairwise distances and then divide by n(n-1)
        let n = align_length as i32;
        let dist_sum = distances.iter().sum::<i32>();
        let np = ((n * (n - 1)) as f32) / 2.0;
        dist_sum as f32 / np
    }
    // 2 given two string sequences, output the hamming distance between them
    // ignore base pair comparisons where there is one or more dash, or N.

    fn hamming_distance(seq1: &String, seq2: &String) -> i32 {
        let mut distance = 0i32;

        for (b1, b2) in seq1.as_bytes().iter().zip(seq2.as_bytes().iter()) {
            if b1 != b2 && *b1 != 45u8 && *b2 != 45u8 && *b1 != 78u8 && *b2 != 78u8 {
                distance += 1;
            } else {
                continue;
            }
        }
        distance
    }

    fn plot_pi_windows(data: Vec<(usize, usize, f32)>) -> Result<(), Box<dyn std::error::Error>> {
        // dimensions of the plot
        let dims = (1280, 2 * 480);
        // x from zero to length vec
        // y from zero to max vec.

        let ymax = data
            .iter()
            .map(|(_a, _b, c)| c)
            .fold(0.0f32, |max, &val| if val > max { val } else { max });
        let xmax = data[data.len() - 1].0;

        let root = BitMapBackend::new("./test.png", (dims.0, dims.1)).into_drawing_area();
        root.fill(&WHITE)?;
        let root = root.margin(10, 10, 10, 10);
        // After this point, we should be able to draw construct a chart context
        let mut chart = ChartBuilder::on(&root)
            .x_label_area_size(20)
            .y_label_area_size(40)
            .set_label_area_size(LabelAreaPosition::Left, (8).percent())
            .set_label_area_size(LabelAreaPosition::Bottom, (8).percent())
            // zero first value needs to be replaced by a minimum I think... for the Y
            .build_cartesian_2d(0f32..xmax as f32, 0f32..ymax as f32)?;

        chart
            .configure_mesh()
            .y_desc("Nucleotide diversity")
            .x_desc("Length along alignment")
            .label_style(TextStyle::from(("sans-serif", 25)))
            .draw()?;

        let line_data: Vec<(f32, f32)> = data.iter().map(|(a, _b, c)| (*a as f32, *c)).collect();

        chart.draw_series(LineSeries::new(line_data, &BLACK))?;

        Ok(())
    }
}
