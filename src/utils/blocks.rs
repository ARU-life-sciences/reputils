use crate::utils::alignment::Alignment;
use std::fmt;

#[derive(Copy, Clone)]
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
        let removed_isolates = Self::remove_isolates(self);
        let start = removed_isolates.0[0].position;
        let end = removed_isolates.0[removed_isolates.0.len() - 1].position;

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

    // this is a little hacky
    // removes elements of a vector if
    // there are no neighbours within 10bp
    // important in removing outliers.

    fn remove_isolates(&self) -> Self {
        let mut index = 0;
        let len = self.0.len();

        let mut out_vec = Vec::new();

        loop {
            if index == len - 1 {
                break;
            }

            let val = *self.0.get(index).unwrap(); // should never panic.
            let left_val;
            if index == 0 {
                left_val = self.0.get(index);
            } else {
                left_val = self.0.get(index - 1);
            }

            let right_val = self.0.get(index + 1);

            // some logic
            match left_val {
                Some(lv) => {
                    match right_val {
                        Some(rv) => {
                            // lv & rv exist
                            if !(val.position - 10 > lv.position)
                                || !(val.position + 10 < rv.position)
                            {
                                out_vec.push(val)
                            }
                        }
                        None => {
                            // lv only exists
                            if !(val.position - 10 > lv.position) {
                                out_vec.push(val)
                            }
                        }
                    }
                }
                None => {
                    match right_val {
                        Some(rv) => {
                            // rv only exists
                            if !(val.position + 10 < rv.position) {
                                out_vec.push(val)
                            }
                        }
                        None => out_vec.push(val),
                    }
                }
            }
            index += 1;
        }
        BlockRecords(out_vec)
    }
}
