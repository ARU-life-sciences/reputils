use crate::utils::alignment::{Alignment, Sequence};
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

    // add internal switch, so we can use the inner API for html
    pub fn trim(
        &self,
        alignment: Alignment,
        extend: usize,
        next_hit: usize,
        internal: bool,
    ) -> Option<Alignment> {
        let removed_isolates = Self::remove_isolates(self, next_hit);
        let start: usize;
        let end: usize;

        // in good cases, the removed isolates vec will be
        match removed_isolates.0.is_empty() {
            true => {
                // start at the beginning
                start = 0;
                // assume the end of the alignment is the end.
                end = alignment.matrix[0].len();
            }
            false => {
                start = removed_isolates.0[0].position;
                end = removed_isolates.0[removed_isolates.0.len() - 1].position;
            }
        }

        // if we want to render a html doc.
        let mut _internal_alignment = Alignment::new();

        for seq in alignment.matrix {
            // more guarding needed here
            let seq_length = seq.sequence.len();
            // calculate new starts and ends
            let new_start: usize;
            let new_end: usize;

            if start - extend <= 0 {
                new_start = 0;
            } else {
                new_start = start - extend;
            }

            if end + extend >= seq_length {
                new_end = seq_length;
            } else {
                new_end = end + extend;
            }

            let trimmed_seq = seq.sequence.get(new_start..new_end);

            if internal {
                // push to internal parseable struct
                match trimmed_seq {
                    Some(t) => _internal_alignment.add_sequence(Sequence {
                        name: seq.name,
                        sequence: t.to_vec(),
                    }),
                    None => {
                        eprintln!("[-]\tHmm, shouldn't ever reach this message. Report as bug!")
                    }
                }
            } else {
                match trimmed_seq {
                    Some(t) => {
                        println!(">{}\n{}", seq.name, std::str::from_utf8(t).unwrap_or(""))
                    }
                    None => {
                        eprintln!("[-]\tHmm, shouldn't ever reach this message. Report as bug!")
                    }
                }
            }
        }
        Some(_internal_alignment)
    }

    // this is a little hacky
    // removes elements of a vector if
    // there are no neighbours within 10bp
    // important in removing outliers.

    // this probably needs to be re-done to exculde any
    // non-consecutive hits.

    fn remove_isolates(&self, next_hit: usize) -> Self {
        let mut index = 0;
        let len = self.0.len();

        let mut out_vec = Vec::new();

        loop {
            if index == len - 1 {
                break;
            }

            let val = *self.0.get(index).unwrap(); // should never panic.
            let left_val;
            // guard against indexing out of range
            if index == 0 {
                left_val = self.0.get(index);
            } else {
                left_val = self.0.get(index - 1);
            }

            // no need to guard here as the break should have it covered
            let right_val = self.0.get(index + 1);

            // some logic
            match left_val {
                Some(lv) => {
                    match right_val {
                        Some(rv) => {
                            // lv & rv exist
                            // if the position of the current hit
                            // minus/plus 10bp is not greater/less (i.e. it includes)
                            // than the previous/next hit, add this hit.
                            if !(val.position - next_hit > lv.position)
                                && !(val.position + next_hit < rv.position)
                            {
                                out_vec.push(val)
                            }
                        }
                        None => {
                            // lv only exists
                            if !(val.position - next_hit > lv.position) {
                                out_vec.push(val)
                            }
                        }
                    }
                }
                None => {
                    match right_val {
                        Some(rv) => {
                            // rv only exists
                            if !(val.position + next_hit < rv.position) {
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
