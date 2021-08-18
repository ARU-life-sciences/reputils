pub mod ttc {
    // trim to consensus
    // trim to core sequence?
    // time to cry?

    use bio::io::fasta;
    use clap::value_t;

    use crate::utils::alignment::{Alignment, Sequence};

    pub fn ttc(matches: &clap::ArgMatches) {
        let fasta = matches.value_of("fasta").unwrap();
        let extend = value_t!(matches.value_of("extend"), usize).unwrap_or_else(|e| e.exit());
        let next_hit = value_t!(matches.value_of("next_hit"), usize).unwrap_or_else(|e| e.exit());
        let miss = value_t!(matches.value_of("missing"), f64).unwrap_or_else(|e| e.exit());
        let iden = value_t!(matches.value_of("identity"), f64).unwrap_or_else(|e| e.exit());

        let reader = fasta::Reader::from_file(fasta).expect("[-]\tPath invalid.");

        // read the fasta into our struct
        let mut matrix = Alignment::new();
        for record in reader.records() {
            let fasta_record = record.expect("[-]\tError during fasta record parsing.");
            matrix.add_sequence(Sequence {
                name: fasta_record.id().to_string(),
                sequence: fasta_record.seq().to_vec(),
            })
        }

        // find the blocks and trim
        let blocks = matrix.find_blocks(miss, iden);
        blocks.trim(matrix, extend, next_hit, false);
    }
}
