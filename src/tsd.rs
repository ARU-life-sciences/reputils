// after an alignment has been trimmed

// TODO: add options for length of each end looked at

pub mod tsd {
    use crate::utils::alignment::{Alignment, Sequence};
    use bio::io::fasta;
    use clap::value_t;

    pub fn find_tsds(matches: &clap::ArgMatches) {
        // for each sequence
        // look at the 20 bases(?) at either end
        // iterate over in windows of length 2 - 12 (2 = 19, 3 = 18, 4 = 17... 12 = 9)
        // store in a hashmap (left and right for each sequence.)
        // 10 hashmap pairs (x no. seq)
        // then ask which hashmap keys are equal
        // of those equal pairs

        let fasta = matches.value_of("fasta").unwrap();
        let length = value_t!(matches.value_of("length"), usize).unwrap_or_else(|e| e.exit());
        let min_window = value_t!(matches.value_of("minimum"), usize).unwrap_or_else(|e| e.exit());
        let max_window = value_t!(matches.value_of("maximum"), usize).unwrap_or_else(|e| e.exit());

        let reader = fasta::Reader::from_file(fasta).expect("[-]\tPath invalid.");

        let mut alignment = Alignment::new();
        for record in reader.records() {
            let record = record.expect("[-]\tError during fasta record parsing.");

            alignment.add_sequence(Sequence {
                name: record.id().to_string(),
                sequence: record.seq().to_vec(),
            })
        }
        // prints a table
        alignment
            .to_tsd_hash(length, min_window, max_window)
            .merge();
    }
}
