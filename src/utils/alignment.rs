/// Module for handling alignments in fasta format
use crate::utils::blocks::{BlockRecord, BlockRecords};
use std::collections::HashMap;
use std::fmt;

/// Sequence is a struct of the raw sequence
/// encoding the sequences as a Vec<u8>
/// and the headers as Strings.

#[derive(Debug)]
pub struct Sequence {
    pub name: String,
    pub sequence: Vec<u8>,
}

impl Sequence {
    pub fn len(&self) -> usize {
        self.sequence.len()
    }
    // an is_empty method for Sequence
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

/// Alignment loads all the Sequences
/// in the fasta file into memory
#[derive(Debug)]
pub struct Alignment {
    pub matrix: Vec<Sequence>,
}

impl Alignment {
    pub fn new() -> Self {
        Alignment { matrix: Vec::new() }
    }
    // add a sequence to the alignment in memory
    pub fn add_sequence(&mut self, seq: Sequence) {
        self.matrix.push(seq)
    }

    pub fn clear(&mut self) {
        self.matrix.clear()
    }

    pub fn find_blocks(&self, miss: f64, iden: f64) -> BlockRecords {
        let t = Self::transpose(self);

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
                .max_by(|a, b| a.1.cmp(b.1))
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
            if per_missing < miss && per_identity > iden {
                blocks.add_record(BlockRecord {
                    position,
                    identity: per_identity,
                    missing: per_missing,
                })
            }
        }
        // println!("{}", blocks);
        blocks
    }

    // transpose the alignment
    fn transpose(&self) -> Vec<Vec<u8>> {
        let v = &self.matrix;
        assert!(!v.is_empty());
        (0..v[0].len())
            .map(|i| v.iter().map(|inner| inner.sequence[i]).collect::<Vec<u8>>())
            .collect()
    }

    pub fn div_windows(
        &self,
        window_size: usize,
        window_step: usize,
        internal: bool,
    ) -> Vec<(usize, usize, f32)> {
        // let mut window = Vec::<String>::new(); // vector to store window
        let mut window = Alignment::new();
        let mut start: usize = 0; // index to start the window
        let mut end: usize = start + window_size; // index to end the window

        let mut data = Vec::new();

        while end <= self.matrix[0].len() {
            for row in &self.matrix {
                let seq = &row.sequence[start..end];
                window.add_sequence(Sequence {
                    name: row.name.clone(),
                    sequence: seq.to_vec(),
                });
            }

            let pi = window.calculate_pi();
            // window contains a vec of chunks.
            if !internal {
                println!("{}\t{}\t{:.3}", start, end, pi);
            }
            data.push((start, end, pi));

            start += window_step;
            end += window_step;

            window.clear();
        }

        data
    }

    // heavily poached from https://github.com/noahaus/sliding-window-scripts/blob/13d872b379a3501ca9b506f8bcccf89a9cd81c8d/tjd/src/main.rs

    pub fn calculate_pi(&self) -> f32 {
        // number of sequences being analyzed
        let align_length = self.matrix.len();
        // vector that holds all the values of the pairwise distances
        let mut distances = Vec::new();

        // calculate hamming distance between two sequences at a time
        for i in 0..align_length {
            for j in 0..align_length {
                if i < j {
                    distances.push(Self::hamming_distance(
                        &self.matrix[i].sequence,
                        &self.matrix[j].sequence,
                    ));
                } else {
                    continue;
                }
            }
        }

        // sum the pairwise distances and then divide by n(n-1)
        let n = align_length as i32;
        let dist_sum = distances.iter().sum::<i32>();
        let np = ((n * (n - 1)) as f32) / 2.0;
        dist_sum as f32 / np
    }

    // 2 given two string sequences, output the hamming distance between them
    // ignore base pair comparisons where there is one or more dash, or N.

    fn hamming_distance(seq1: &[u8], seq2: &[u8]) -> i32 {
        let mut distance = 0i32;

        for (b1, b2) in seq1.iter().zip(seq2.iter()) {
            if b1 != b2 && *b1 != 45u8 && *b2 != 45u8 && *b1 != 78u8 && *b2 != 78u8 {
                distance += 1;
            } else {
                continue;
            }
        }
        distance
    }

    // Get frequencies of each nucleotide at each position in an alignment (profile)

    pub fn get_profile(&self) -> Vec<HashMap<u8, usize>> {
        // as all sequence lengths should be the same
        let sequence_length = self.matrix[0].sequence.len();
        let mut profile = Vec::with_capacity(sequence_length);
        for i in 0..sequence_length {
            let position_string = self
                .matrix
                .iter()
                .map(|record| *record.sequence.get(i).unwrap())
                .collect::<Vec<u8>>();
            let mut profile_hash = HashMap::new();
            for base in position_string {
                *profile_hash.entry(base).or_insert(0) += 1;
            }
            profile.push(profile_hash);
        }
        profile
    }

    pub fn to_tsd_hash(&self, length: usize, min_window: usize, max_window: usize) -> TSDHash {
        let sequence_length = self.matrix[0].sequence.len();

        let mut tsd_hash = TSDHash::new();

        for record in &self.matrix {
            // FIX THESE UNWRAPS
            let left_seq = record.sequence.get(0..length).unwrap();
            let right_seq = record
                .sequence
                .get(sequence_length - length..sequence_length - 1)
                .unwrap();

            for window in min_window..=max_window {
                let l_window_iterator = left_seq.windows(window);
                let r_window_iterator = right_seq.windows(window);

                let mut left = HashMap::new();
                let mut right = HashMap::new();

                for l_w in l_window_iterator {
                    if !l_w.contains(&45) {
                        *left.entry(l_w.to_vec()).or_insert(0) += 1usize;
                    }
                }
                for r_w in r_window_iterator {
                    if !r_w.contains(&45) {
                        *right.entry(r_w.to_vec()).or_insert(0) += 1usize;
                    }
                }

                tsd_hash.add_entry(AlignHash {
                    name: record.name.clone(),
                    window,
                    left,
                    right,
                })
            }
        }
        tsd_hash
    }
}

// query alignment ends
pub struct AlignHash {
    pub name: String,
    pub window: usize,
    pub left: HashMap<Vec<u8>, usize>,
    pub right: HashMap<Vec<u8>, usize>,
}

pub struct TSDHash(Vec<AlignHash>);

impl TSDHash {
    pub fn new() -> Self {
        TSDHash(Vec::new())
    }

    pub fn add_entry(&mut self, other: AlignHash) {
        self.0.push(AlignHash {
            name: other.name,
            window: other.window,
            left: other.left,
            right: other.right,
        })
    }

    // some more work to do here;
    // if the smaller 'tsd's' are substrings of a larger
    // they can be removed (they're redundant.)
    pub fn merge(&self, internal: bool) -> Option<HashMap<String, Vec<String>>> {
        let mut res = Vec::new();
        for hash in &self.0 {
            // merge left and right hashmaps
            let merged = Self::merge_hashmaps(&hash.left, &hash.right);

            res.push((hash.name.clone(), merged));
        }

        // group by ID
        let mut grouped_id = HashMap::new();
        for (id, map) in res {
            // collect into a vec
            let tsds = map.keys().collect::<Vec<&Vec<u8>>>();
            // turn to strings
            let mut tsds_str = tsds
                .iter()
                .map(|e| std::str::from_utf8(e).unwrap().to_owned())
                .collect::<Vec<String>>();

            grouped_id
                .entry(id)
                .or_insert(Vec::new())
                .append(&mut tsds_str);
        }

        // collect into another hashmap to remove substrings
        // long live immutability!
        let mut grouped_id_no_substr = HashMap::new();

        for (k, v) in grouped_id {
            let mut no_substr = Vec::new();

            for s1_i in 0..v.len() {
                let mut substring_found = false;
                for s2 in &v {
                    if &v[s1_i] != s2 && s2.contains(&v[s1_i]) {
                        substring_found = true;
                    }
                }
                if !substring_found {
                    no_substr.push(v[s1_i].clone());
                }
            }
            grouped_id_no_substr.insert(k, no_substr);
        }

        if internal {
            Some(grouped_id_no_substr)
        } else {
            // for command line printing only.
            println!("Table of potential TSD's");
            println!("ID\tLength\tTSD's");
            for (id, tsds) in grouped_id_no_substr {
                if !tsds.is_empty() {
                    println!("{}\t{}", id, tsds.join("\t"));
                }
            }
            None
        }
    }

    // take TSDHash,
    fn merge_hashmaps(
        map1: &HashMap<Vec<u8>, usize>,
        map2: &HashMap<Vec<u8>, usize>,
    ) -> HashMap<Vec<u8>, usize> {
        let mut merged = HashMap::new();

        for (k1, v1) in map1 {
            for (k2, v2) in map2 {
                if k1 == k2 {
                    merged.insert(k1.clone(), v1 + v2);
                }
            }
        }
        merged
    }
}

impl fmt::Display for TSDHash {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for record in &self.0 {
            writeln!(
                f,
                "Sequence: {}\nWindow size: {}\nLeft hash: {:?}\nRight hash: {:?}\n",
                record.name, record.window, record.left, record.right
            )?;
        }
        Ok(())
    }
}
