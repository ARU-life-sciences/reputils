pub mod div {
    use bio::io::fasta;
    use clap::value_t;
    use plotters::prelude::*;

    use crate::utils::alignment::{Alignment, Sequence};

    pub fn diversity_windows(matches: &clap::ArgMatches) {
        let fasta = matches.value_of("fasta").unwrap();
        let window_size = value_t!(matches.value_of("window"), usize).unwrap_or_else(|e| e.exit());
        let window_step = value_t!(matches.value_of("step"), usize).unwrap_or_else(|e| e.exit());
        let plot = matches.is_present("plot");
        let dir = value_t!(matches.value_of("dir"), String).unwrap_or_else(|e| e.exit());
        let name = value_t!(matches.value_of("name"), String).unwrap_or_else(|e| e.exit());

        let mut reader = fasta::Reader::from_file(fasta)
            .expect("[-]\tPath invalid.")
            .records();

        let mut alignment = Alignment::new();

        while let Some(Ok(record)) = reader.next() {
            // alignment.push(std::str::from_utf8(record.seq()).unwrap().to_string());
            alignment.add_sequence(Sequence {
                name: record.id().to_string(),
                sequence: record.seq().to_vec(),
            });
        }

        // let mut window = Vec::<String>::new(); // vector to store window
        let mut window = Alignment::new();
        let mut start: usize = 0; // index to start the window
        let mut end: usize = start + window_size; // index to end the window

        let mut data = Vec::new();

        while end <= alignment.matrix[0].len() {
            for row in &alignment.matrix {
                let seq = &row.sequence[start..end];
                window.add_sequence(Sequence {
                    name: row.name.clone(),
                    sequence: seq.to_vec(),
                });
            }

            let pi = window.calculate_pi();
            // window contains a vec of chunks.
            println!("{}\t{}\t{:.3}", start, end, pi);
            data.push((start, end, pi));

            start += window_step;
            end += window_step;

            window.matrix.clear();
        }

        if plot {
            // do plot
            plot_pi_windows(data, &dir, &name).expect("Couldn't make the plot :(");
        }
    }

    fn plot_pi_windows(
        data: Vec<(usize, usize, f32)>,
        dir: &str,
        name: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        // dimensions of the plot
        let dims = (1280, 2 * 480);
        // x from zero to length vec
        // y from zero to max vec.

        let ymax = data
            .iter()
            .map(|(_a, _b, c)| c)
            .fold(0.0f32, |max, &val| if val > max { val } else { max });
        let xmax = data[data.len() - 1].0;

        let path = format!("{}/{}.png", dir, name);

        let root = BitMapBackend::new(&path, (dims.0, dims.1)).into_drawing_area();
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
