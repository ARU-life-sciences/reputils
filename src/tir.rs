pub mod tir {

    use bio::alignment::pairwise::*;
    use bio::io::fasta;

    use crate::con::con::get_consensus;
    use crate::utils::alignment::{Alignment, Sequence};
    use crate::utils::revcomp::revcomp::reverse_complement;
    use crate::utils::seqcount::SequentialCount;

    // use this module to perform an alignment between consensus and itself (but revcomp)
    // might be swift evidence of a TIR.
    // TODO: might be useful to trim the pretty alignment? Here or in the HTML.

    // maybe this should take the whole alignment, trim (optional), make consensus, then do self alignment.

    pub fn revcomp_alignment(matches: &clap::ArgMatches) {
        let fasta = matches.value_of("fasta").unwrap();
        let show_alignment = matches.is_present("show");

        // read in the fasta from file
        let mut reader = fasta::Reader::from_file(fasta)
            .expect("[-]\tPath invalid.")
            .records();

        let mut read_number = 0i32;

        let mut alignment = Alignment::new();

        while let Some(Ok(record)) = reader.next() {
            // alignment.push(std::str::from_utf8(record.seq()).unwrap().to_string());
            alignment.add_sequence(Sequence {
                name: record.id().to_string(),
                sequence: record.seq().to_vec(),
            });
            read_number += 1;
        }

        let alignment_blocks = alignment.find_blocks(0.1, 0.8);
        let alignment = alignment_blocks.trim(alignment, 15, 1, true);

        let profile = alignment.unwrap().get_profile();
        let mut forward_consensus = get_consensus(profile, read_number);

        // remove all gaps as they mess up the alignment
        forward_consensus.retain(|&e| e != 45);

        let reverse_consensus = reverse_complement(&forward_consensus);

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        // gap open score: -5, gap extension score: -1
        let mut aligner = Aligner::with_capacity(
            forward_consensus.len(),
            reverse_consensus.len(),
            -5,
            -1,
            &score,
        );
        let alignment = aligner.semiglobal(&forward_consensus, &reverse_consensus);
        let alignment_pretty = alignment.pretty(&forward_consensus, &reverse_consensus);

        let mut lim = 0;
        let mut start = 0;
        let mut end = 0;

        for (base, count) in SequentialCount::new(alignment.operations.iter()) {
            end += count;
            if lim < 10 {
                println!("{} - {}: {:?} occurs {} times", start, end, base, count);
            } else {
                break;
            }
            lim += 1;
            start = end;
        }

        if show_alignment {
            print!("\n{}", alignment_pretty);
        }
    }
}
