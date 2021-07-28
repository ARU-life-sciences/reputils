pub mod tir {

    use bio::alignment::pairwise::*;
    use bio::io::fasta;
    use std::iter::Peekable;
    // use this module to perform an alignment between consensus and itself (but revcomp)
    // might be swift evidence of a TIR.

    pub fn revcomp_alignment(matches: &clap::ArgMatches) {
        let fasta = matches.value_of("fasta").unwrap();
        let show_alignment = matches.is_present("show");
        // read number check
        let mut read_number = 0i32;
        let mut forward: Vec<u8> = Vec::new();
        // read in the fasta from file
        let mut reader = fasta::Reader::from_file(fasta)
            .expect("[-]\tPath invalid.")
            .records();
        while let Some(Ok(record)) = reader.next() {
            forward = record.seq().to_owned();
            read_number += 1;
        }

        if read_number > 1 {
            panic!("[-]\tThere should only be one sequence in the fasta file.")
        }

        let reverse = reverse_complement(&forward);

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        // gap open score: -5, gap extension score: -1
        let mut aligner = Aligner::with_capacity(forward.len(), reverse.len(), -5, -1, &score);
        let alignment = aligner.semiglobal(&forward, &reverse);
        let alignment_pretty = alignment.pretty(&forward, &reverse);

        let mut lim = 0;
        let mut start = 0;
        let mut end = 0;

        for (base, count) in SequentialCount::new(alignment.operations.iter()) {
            end += count;
            if lim < 10 {
                println!("{} - {}: {:?} occurs {} times", start, end, base, count);
            }
            lim += 1;
            start = end;
        }

        if show_alignment {
            print!("\n{}", alignment_pretty);
        }
    }

    fn reverse_complement(dna: &[u8]) -> Vec<u8> {
        let dna_vec = dna.to_vec();
        let mut revcomp = Vec::new();

        for base in dna_vec.iter() {
            revcomp.push(switch_base(*base))
        }
        revcomp.as_mut_slice().reverse();
        revcomp
    }

    fn switch_base(c: u8) -> u8 {
        match c {
            b'A' => b'T',
            b'a' => b't',
            b'C' => b'G',
            b'c' => b'g',
            b'T' => b'A',
            b't' => b'a',
            b'G' => b'C',
            b'g' => b'c',
            b'N' => b'N',
            b'n' => b'n',
            _ => b'N',
        }
    }

    struct SequentialCount<I>
    where
        I: Iterator,
    {
        iter: Peekable<I>,
    }

    impl<I> SequentialCount<I>
    where
        I: Iterator,
    {
        fn new(iter: I) -> Self {
            SequentialCount {
                iter: iter.peekable(),
            }
        }
    }

    impl<I> Iterator for SequentialCount<I>
    where
        I: Iterator,
        I::Item: Eq,
    {
        type Item = (I::Item, usize);

        fn next(&mut self) -> Option<Self::Item> {
            // Check the next value in the inner iterator
            match self.iter.next() {
                // There is a value, so keep it
                Some(head) => {
                    // We've seen one value so far
                    let mut count = 1;
                    // Check to see what the next value is without
                    // actually advancing the inner iterator
                    while self.iter.peek() == Some(&head) {
                        // It's the same value, so go ahead and consume it
                        self.iter.next();
                        count += 1;
                    }
                    // The next element doesn't match the current value
                    // complete this iteration
                    Some((head, count))
                }
                // The inner iterator is complete, so we are also complete
                None => None,
            }
        }
    }
}
