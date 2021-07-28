use std::collections::HashMap;

pub mod ttc {
    use bio::io::fasta;
    use clap::value_t;

    use super::{Alignment, Sequence};

    pub fn ttc(matches: &clap::ArgMatches) {
        let fasta = matches.value_of("fasta").unwrap();
        let extend = value_t!(matches.value_of("extend"), usize).unwrap_or_else(|e| e.exit());

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
        let blocks = matrix.find_blocks();
        blocks.trim(matrix, extend);
    }
}

/// Sequence is s struct of the raw sequence
/// encoding the sequences as a Vec<u8>
#[derive(Debug)]
pub struct Sequence {
    name: String,
    sequence: Vec<u8>,
}

impl Sequence {
    pub fn len(&self) -> usize {
        self.sequence.len()
    }
}

/// Alignment merges all the Sequences
/// in the fasta file
#[derive(Debug)]
pub struct Alignment {
    matrix: Vec<Sequence>,
}

pub struct BlockRecord {
    pub position: usize,
    pub identity: f64,
    pub missing: f64,
}

pub struct BlockRecords(Vec<BlockRecord>);

use std::fmt;

impl fmt::Display for BlockRecords {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for record in &self.0 {
            writeln!(
                f,
                "Position: {} -- %ID: {} -- %Miss: {}",
                record.position + 1, // add one due to zero indexing.
                record.identity,
                record.missing
            )?;
        }
        Ok(())
    }
}

impl BlockRecords {
    pub fn add_record(&mut self, value: BlockRecord) {
        self.0.push(value)
    }
    pub fn trim(&self, alignment: Alignment, extend: usize) {
        let start = self.0[0].position;
        let end = self.0[self.0.len() - 1].position;

        for seq in alignment.matrix {
            let trimmed_seq = seq.sequence.get(start - extend..end + extend);
            match trimmed_seq {
                Some(t) => {
                    println!(">{}\n{}", seq.name, std::str::from_utf8(t).unwrap_or(""))
                }
                None => {
                    eprintln!("[-]\tRange out of bounds; use smaller value for *extend*.")
                }
            }
        }
    }
}

impl Alignment {
    pub fn new() -> Alignment {
        Alignment { matrix: Vec::new() }
    }
    // add a sequence to the alignment in memory
    pub fn add_sequence(&mut self, seq: Sequence) {
        self.matrix.push(seq)
    }

    pub fn find_blocks(&self) -> BlockRecords {
        let t = Alignment::transpose(&self);

        let mut blocks: BlockRecords = BlockRecords(Vec::new());
        for (position, column) in t.iter().enumerate() {
            let col_len = column.len();
            // calculate identity & missing
            let mut table = HashMap::new();
            for base in column {
                *table.entry(base).or_insert(0) += 1;
            }
            // percent identical. bad error handling here.
            let identity = table
                .iter()
                .max_by(|a, b| a.1.cmp(&b.1))
                .map(|(k, v)| (k, v))
                .unwrap_or((&&0u8, &0i32));
            let per_identity = {
                if identity.0 == &&45u8 {
                    0.0
                } else {
                    *identity.1 as f64 / col_len as f64
                }
            };
            let per_missing = *table.get(&45u8).unwrap_or(&0) as f64 / col_len as f64;

            // add a block IF
            // missing (dashes) below 0.1 AND column should be 0.9 identical.
            if per_missing < 0.1 && per_identity > 0.9 {
                blocks.add_record(BlockRecord {
                    position,
                    identity: per_identity,
                    missing: per_missing,
                })
            }
        }
        blocks
    }

    // transpose the alignment
    fn transpose(&self) -> Vec<Vec<u8>> {
        let v = &self.matrix;
        assert!(!v.is_empty());
        (0..v[0].len())
            .map(|i| {
                v.iter()
                    .map(|inner| inner.sequence[i].clone())
                    .collect::<Vec<u8>>()
            })
            .collect()
    }
}
