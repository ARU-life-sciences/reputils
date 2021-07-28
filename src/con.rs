pub mod con {

    use bio::io::fasta;
    use std::collections::HashMap;
    use std::fmt::{Display, Error, Formatter};

    // with much help from https://github.com/Ninjani/rosalind/blob/master/s_cons/src/lib.rs

    pub fn make_consensus(matches: &clap::ArgMatches) {
        // parse command line args
        let fasta = matches.value_of("fasta").unwrap();
        let name = matches.value_of("name").unwrap();
        let append = matches.is_present("append");

        // do some read length checks.
        // get the number of reads for free.
        let mut check_read_lengths = Vec::new();
        let mut read_number = 0i32;
        // read in the fasta from file
        let mut reader = fasta::Reader::from_file(fasta)
            .expect("[-]\tPath invalid.")
            .records();
        while let Some(Ok(record)) = reader.next() {
            read_number += 1;
            check_read_lengths.push(record.seq().len());
        }
        let check_lens = is_all_same(&check_read_lengths);

        if !check_lens {
            panic!("[-]\tAll sequences in the fasta file are not the same length.")
        }

        let reader = fasta::Reader::from_file(fasta).expect("[-]\tPath invalid.");

        // collect sequences into memory
        // should be fine for small(ish) alignments
        let mut sequences = Vec::new();
        for result in reader.records() {
            let record = result.expect("[-]\tError during fasta record parsing.");
            let id = record.id().to_owned();
            let sequence = record.seq().to_owned();
            sequences.push((id, sequence));
        }

        // containing the frequencies of each nucleotide at each column
        let profile = get_profile(sequences.clone());
        let consensus = get_consensus(profile, read_number);

        if append {
            println!(">{}\n{}", name, WriteSequence(consensus));
        } else {
            for sequence in sequences {
                println!(">{}\n{}", sequence.0, WriteSequence(sequence.1));
            }
            println!(">{}\n{}", name, WriteSequence(consensus));
        }
    }

    // Get frequencies of each nucleotide at each position in a collection of sequences (profile)

    fn get_profile(sequences: Vec<(String, Vec<u8>)>) -> Vec<HashMap<u8, usize>> {
        // as all sequence lengths should be the same
        let sequence_length = sequences[0].1.len();
        let mut profile = Vec::with_capacity(sequence_length);
        for i in 0..sequence_length {
            let position_string = sequences
                .iter()
                .map(|(_id, sequence)| *sequence.iter().nth(i).unwrap())
                .collect::<Vec<u8>>();
            let mut profile_hash = HashMap::new();
            for base in position_string {
                *profile_hash.entry(base).or_insert(0) += 1;
            }
            profile.push(profile_hash);
        }
        profile
    }

    /// Get consensus sequence from a profile
    fn get_consensus(profile: Vec<HashMap<u8, usize>>, read_number: i32) -> Vec<u8> {
        // initiate consensus
        let mut consensus = Vec::with_capacity(profile.len());

        let minority = 0.4;
        // iterate over the vector of hashmaps
        for counts in profile.iter() {
            // frequencies of A, C, G, T's
            let g_c = counts.get(&71).unwrap_or(&0) + counts.get(&103).unwrap_or(&0); // 71 == G; 103 == g
            let c_c = counts.get(&67).unwrap_or(&0) + counts.get(&99).unwrap_or(&0); // 67 == C; 99 == c
            let a_c = counts.get(&65).unwrap_or(&0) + counts.get(&97).unwrap_or(&0); // 65 == A; 97 == a
            let t_c = counts.get(&84).unwrap_or(&0) + counts.get(&116).unwrap_or(&0); // 84 == T; 116 == t

            // deal with single best variants first
            // this is a very crazy way to do this. But seems to work.
            if g_c > c_c && g_c > a_c && g_c > t_c && g_c as f32 > (minority * read_number as f32) {
                consensus.push(71u8);
            } else if c_c > g_c
                && c_c > a_c
                && c_c > t_c
                && c_c as f32 > (minority * read_number as f32)
            {
                consensus.push(67u8);
            } else if a_c > g_c
                && a_c > c_c
                && a_c > t_c
                && a_c as f32 > (minority * read_number as f32)
            {
                consensus.push(65u8);
            } else if t_c > g_c
                && t_c > a_c
                && t_c > c_c
                && t_c as f32 > (minority * read_number as f32)
            {
                consensus.push(84u8);
            }
            // now deal with duals (IUPAC)
            // g/t == k
            else if g_c > c_c
                && g_c > a_c
                && g_c == t_c
                && g_c as f32 > (minority * read_number as f32)
            {
                consensus.push(75u8);
            }
            // g/c == s
            else if g_c > a_c
                && g_c > t_c
                && g_c == c_c
                && g_c as f32 > (minority * read_number as f32)
            {
                consensus.push(83u8);
            }
            // g/a == r
            else if g_c > c_c
                && g_c > t_c
                && g_c == a_c
                && g_c as f32 > (minority * read_number as f32)
            {
                consensus.push(82u8);
            }
            // a/c == m
            else if a_c > g_c
                && a_c > t_c
                && a_c == c_c
                && a_c as f32 > (minority * read_number as f32)
            {
                consensus.push(77u8);
            }
            // a/t == w
            else if a_c > g_c
                && a_c > c_c
                && a_c == t_c
                && a_c as f32 > (minority * read_number as f32)
            {
                consensus.push(87u8);
            }
            // c/t == y
            else if c_c > g_c
                && c_c > a_c
                && c_c == t_c
                && c_c as f32 > (minority * read_number as f32)
            {
                consensus.push(89u8);
            }
            // a/c/g == v
            else if a_c > t_c
                && a_c == c_c
                && a_c == g_c
                && a_c as f32 > (minority * read_number as f32)
            {
                consensus.push(86u8);
            }
            // a/c/t == h
            else if a_c > g_c
                && a_c == c_c
                && a_c == t_c
                && a_c as f32 > (minority * read_number as f32)
            {
                consensus.push(72u8);
            }
            // a/g/t == d
            else if a_c > c_c
                && a_c == g_c
                && a_c == t_c
                && a_c as f32 > (minority * read_number as f32)
            {
                consensus.push(68u8);
            }
            // c/g/t == b
            else if c_c > a_c
                && c_c == g_c
                && c_c == t_c
                && c_c as f32 > (minority * read_number as f32)
            {
                consensus.push(66u8);
            }
            // a/c/g/t == N
            // should this be N????
            else if c_c == a_c
                && c_c == g_c
                && c_c == t_c
                && c_c as f32 > (minority * read_number as f32)
            {
                consensus.push(78u8);
            }
            // add ? if any majority nucleotide frequency below 0.5
            // but if < 0.2, then put a dash
            // these numbers can of course be changed.
            else if (c_c as f32) < (minority * read_number as f32)
                && (c_c as f32) > (0.3 * read_number as f32)
                || (g_c as f32) < (minority * read_number as f32)
                    && (g_c as f32) > (0.3 * read_number as f32)
                || (a_c as f32) < (minority * read_number as f32)
                    && (a_c as f32) > (0.3 * read_number as f32)
                || (t_c as f32) < (minority * read_number as f32)
                    && (t_c as f32) > (0.3 * read_number as f32)
            {
                consensus.push(63u8);
            } else {
                consensus.push(45u8);
            }
        }
        consensus
    }

    fn is_all_same(arr: &[usize]) -> bool {
        if arr.is_empty() {
            return true;
        }
        let first = arr[0];
        arr.iter().all(|&item| item == first)
    }

    // display for Vec<i32>
    #[derive(Clone)]
    pub struct WriteSequence(pub Vec<u8>);

    impl Display for WriteSequence {
        fn fmt(&self, f: &mut Formatter) -> Result<(), Error> {
            let sequence = std::str::from_utf8(&self.0).unwrap();
            write!(f, "{}", sequence)
        }
    }
}
