use crate::utils::alignment::Alignment;
use std::fmt;

pub struct BlockRecord {
    pub position: usize,
    pub identity: f64,
    pub missing: f64,
}

pub struct BlockRecords(pub Vec<BlockRecord>);

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
