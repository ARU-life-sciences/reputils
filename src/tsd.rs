// NOT YET FUNCTIONAL. TODO.

pub mod tsd {
    use bio::io::fasta;
    use bio::pattern_matching::shift_and;
    use std::collections::HashMap;
    use std::ops::{Add, Div, Sub};
    extern crate stats;

    pub fn find_tsds(matches: &clap::ArgMatches) {
        // for each sequence
        // isolate the sequences around aligned chunks
        // split into chunks
        // if chunk contains 2 or more dashes, N's or ?'s, remove
        // if chunks unique, remove
        // HashMap<String, Vec<i32>>, where i32 is the position
        // if position in the vec ^ is within 300bp, remove
        // use information from other sequences to determine TSD (i.e. those with >= 2 entries in the hashmap)
        // then print out possible sequences

        let fasta = matches.value_of("fasta").unwrap();

        // length of the TSD
        let tsd_length = 3;
        // length either side of a match (possible start/end of alignment)
        let length = 3;

        // get index of a dash and a nucelotide
        // extract 20 nucleotides around this point
        let dash_a = b"-A";
        let dash_g = b"-G";
        let dash_c = b"-C";
        let dash_t = b"-T";

        let a_dash = b"A-";
        let g_dash = b"G-";
        let c_dash = b"C-";
        let t_dash = b"T-";

        let shiftand_a_l = shift_and::ShiftAnd::new(dash_a);
        let shiftand_g_l = shift_and::ShiftAnd::new(dash_g);
        let shiftand_c_l = shift_and::ShiftAnd::new(dash_c);
        let shiftand_t_l = shift_and::ShiftAnd::new(dash_t);

        let shiftand_a_r = shift_and::ShiftAnd::new(a_dash);
        let shiftand_g_r = shift_and::ShiftAnd::new(g_dash);
        let shiftand_c_r = shift_and::ShiftAnd::new(c_dash);
        let shiftand_t_r = shift_and::ShiftAnd::new(t_dash);

        let reader_search = fasta::Reader::from_file(fasta).expect("[-]\tPath invalid.");

        // vec to collect all of the sequences
        let mut substrs = Vec::new();

        for record in reader_search.records() {
            let record = record.expect("[-]\tError during fasta record parsing.");
            let sequence = record.seq().to_ascii_uppercase();

            let occ_a_l = shiftand_a_l.find_all(&sequence);
            let occ_g_l = shiftand_g_l.find_all(&sequence);
            let occ_c_l = shiftand_c_l.find_all(&sequence);
            let occ_t_l = shiftand_t_l.find_all(&sequence);

            let occ_a_r = shiftand_a_r.find_all(&sequence);
            let occ_g_r = shiftand_g_r.find_all(&sequence);
            let occ_c_r = shiftand_c_r.find_all(&sequence);
            let occ_t_r = shiftand_t_r.find_all(&sequence);

            for i in occ_a_l {
                let mut index1 = i;
                let mut index2 = i;
                if i < length {
                    index1 = length;
                } else if i > &sequence.len() - length {
                    index2 = &sequence.len() - length;
                }
                let substr = &record.seq()[index1 - length..index2 + length];
                substrs.push((i, std::str::from_utf8(substr).unwrap().to_owned()));
            }
            for i in occ_g_l {
                let mut index1 = i;
                let mut index2 = i;
                if i < length {
                    index1 = length;
                } else if i > &sequence.len() - length {
                    index2 = &sequence.len() - length;
                }
                let substr = &record.seq()[index1 - length..index2 + length];
                substrs.push((i, std::str::from_utf8(substr).unwrap().to_owned()));
            }
            for i in occ_c_l {
                let mut index1 = i;
                let mut index2 = i;
                if i < length {
                    index1 = length;
                } else if i > &sequence.len() - length {
                    index2 = &sequence.len() - length;
                }
                let substr = &record.seq()[index1 - length..index2 + length];
                substrs.push((i, std::str::from_utf8(substr).unwrap().to_owned()));
            }
            for i in occ_t_l {
                let mut index1 = i;
                let mut index2 = i;
                if i < length {
                    index1 = length;
                } else if i > &sequence.len() - length {
                    index2 = &sequence.len() - length;
                }
                let substr = &record.seq()[index1 - length..index2 + length];
                substrs.push((i, std::str::from_utf8(substr).unwrap().to_owned()));
            }
            for i in occ_a_r {
                let mut index1 = i;
                let mut index2 = i;
                if i < length {
                    index1 = length;
                } else if i > &sequence.len() - length {
                    index2 = &sequence.len() - length;
                }
                let substr = &record.seq()[index1 - length..index2 + length];
                substrs.push((i, std::str::from_utf8(substr).unwrap().to_owned()));
            }
            for i in occ_g_r {
                let mut index1 = i;
                let mut index2 = i;
                if i < length {
                    index1 = length;
                } else if i > &sequence.len() - length {
                    index2 = &sequence.len() - length;
                }
                let substr = &record.seq()[index1 - length..index2 + length];
                substrs.push((i, std::str::from_utf8(substr).unwrap().to_owned()));
            }
            for i in occ_c_r {
                let mut index1 = i;
                let mut index2 = i;
                if i < length {
                    index1 = length;
                } else if i > &sequence.len() - length {
                    index2 = &sequence.len() - length;
                }
                let substr = &record.seq()[index1 - length..index2 + length];
                substrs.push((i, std::str::from_utf8(substr).unwrap().to_owned()));
            }
            for i in occ_t_r {
                let mut index1 = i;
                let mut index2 = i;
                if i < length {
                    index1 = length;
                } else if i > &sequence.len() - length {
                    index2 = &sequence.len() - length;
                }
                let substr = &record.seq()[index1 - length..index2 + length];
                substrs.push((i, std::str::from_utf8(substr).unwrap().to_owned()));
            }
        }

        let mut window_hash: HashMap<String, Vec<usize>> = HashMap::new();
        for (i, substr) in substrs {
            let str_substr = substr.as_str().as_bytes();
            let windows = str_substr.windows(tsd_length);

            for window in windows {
                if window.contains(&45) {
                    continue;
                }
                window_hash
                    .entry(std::str::from_utf8(window).unwrap().to_owned())
                    .or_insert(Vec::new())
                    .push(i);
            }
        }
        let vals: Vec<Vec<usize>> = window_hash.values().cloned().collect();
        let vals_flat = flatten(vals);

        // and make a histogram
        // want breaks to be around 50bp?
        let mut hist = Histogram::with_buckets(250);

        for i in vals_flat {
            hist.add(i as u64);
        }

        println!("{}", hist);

        // there might be something in this.
        // in the distribution, find the top two regions with the lowest matches
        //
    }

    fn flatten<T>(nested: Vec<Vec<T>>) -> Vec<T> {
        nested.into_iter().flatten().collect()
    }

    use std::cmp;
    use std::collections::btree_map::Range;
    use std::collections::BTreeMap;
    use std::fmt;

    /// A histogram is a collection of samples, sorted into buckets.
    ///
    /// See the crate level documentation for more details.
    #[derive(Debug, Clone)]
    pub struct Histogram {
        num_buckets: u64,
        samples: BTreeMap<u64, u64>,
        stats: stats::OnlineStats,
        minmax: stats::MinMax<u64>,
    }

    impl Histogram {
        /// Construct a new histogram with the given number of buckets.
        ///
        /// ## Panics
        ///
        /// Panics if the number of buckets is zero.
        pub fn with_buckets(num_buckets: u64) -> Histogram {
            assert!(num_buckets > 0);
            Histogram {
                num_buckets,
                samples: Default::default(),
                stats: Default::default(),
                minmax: Default::default(),
            }
        }

        /// Add a new sample to this histogram.
        pub fn add(&mut self, sample: u64) {
            *self.samples.entry(sample).or_insert(0) += 1;
            self.minmax.add(sample);
            self.stats.add(sample);
        }

        /// Get an iterator over this histogram's buckets.
        pub fn buckets(&self) -> Buckets {
            Buckets {
                histogram: self,
                index: 0,
            }
        }
    }

    impl fmt::Display for Histogram {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            use std::fmt::Write;

            let num_samples: u64 = self.samples.values().sum();
            writeln!(f, "# Number of samples = {}", num_samples)?;
            if num_samples == 0 {
                return Ok(());
            }

            let min = self.minmax.min().unwrap();
            let max = self.minmax.max().unwrap();

            writeln!(f, "# Min = {}", min)?;
            writeln!(f, "# Max = {}", max)?;
            writeln!(f, "#")?;

            let mean = self.stats.mean();
            let dev = self.stats.stddev();
            let var = self.stats.variance();

            writeln!(f, "# Mean = {}", mean)?;
            writeln!(f, "# Standard deviation = {}", dev)?;
            writeln!(f, "# Variance = {}", var)?;
            writeln!(f, "#")?;

            let max_bucket_count = self.buckets().map(|b| b.count()).fold(0, cmp::max);

            const WIDTH: u64 = 50;
            let count_per_char = cmp::max(max_bucket_count / WIDTH, 1);

            writeln!(f, "# Each ∎ is a count of {}", count_per_char)?;
            writeln!(f, "#")?;

            let mut count_str = String::new();

            let widest_count = self.buckets().fold(0, |n, b| {
                count_str.clear();
                write!(&mut count_str, "{}", b.count()).unwrap();
                cmp::max(n, count_str.len())
            });

            let mut end_str = String::new();
            let widest_range = self.buckets().fold(0, |n, b| {
                end_str.clear();
                write!(&mut end_str, "{}", b.end()).unwrap();
                cmp::max(n, end_str.len())
            });

            let mut start_str = String::with_capacity(widest_range);

            for bucket in self.buckets() {
                start_str.clear();
                write!(&mut start_str, "{}", bucket.start()).unwrap();
                for _ in 0..widest_range - start_str.len() {
                    start_str.insert(0, ' ');
                }

                end_str.clear();
                write!(&mut end_str, "{}", bucket.end()).unwrap();
                for _ in 0..widest_range - end_str.len() {
                    end_str.insert(0, ' ');
                }

                count_str.clear();
                write!(&mut count_str, "{}", bucket.count()).unwrap();
                for _ in 0..widest_count - count_str.len() {
                    count_str.insert(0, ' ');
                }

                write!(f, "{} .. {} [ {} ]: ", start_str, end_str, count_str)?;
                for _ in 0..bucket.count() / count_per_char {
                    write!(f, "∎")?;
                }
                writeln!(f)?;
            }

            Ok(())
        }
    }

    /// An iterator over the buckets in a histogram.
    #[derive(Debug, Clone)]
    pub struct Buckets<'a> {
        histogram: &'a Histogram,
        index: u64,
    }

    impl<'a> Iterator for Buckets<'a> {
        type Item = Bucket<'a>;

        fn next(&mut self) -> Option<Self::Item> {
            if self.index >= self.histogram.num_buckets {
                return None;
            }

            let (min, max) = match (self.histogram.minmax.min(), self.histogram.minmax.max()) {
                (Some(&min), Some(&max)) => (min, max),
                _ => return None,
            };

            let range = max - min;
            let range = range + (range % self.histogram.num_buckets);

            let bucket_size = range / self.histogram.num_buckets;
            let bucket_size = cmp::max(1, bucket_size);

            let start = min + self.index * bucket_size;
            let end = min + (self.index + 1) * bucket_size;

            self.index += 1;

            Some(Bucket {
                start,
                end,
                range: if self.index == self.histogram.num_buckets {
                    self.histogram.samples.range(start..)
                } else {
                    self.histogram.samples.range(start..end)
                },
            })
        }
    }

    /// A bucket is a range of samples and their count.
    #[derive(Clone)]
    pub struct Bucket<'a> {
        start: u64,
        end: u64,
        range: Range<'a, u64, u64>,
    }

    impl<'a> Bucket<'a> {
        /// The number of samples in this bucket's range.
        pub fn count(&self) -> u64 {
            self.range.clone().map(|(_, count)| count).sum()
        }

        /// The start of this bucket's range.
        pub fn start(&self) -> u64 {
            self.start
        }

        /// The end of this bucket's range.
        pub fn end(&self) -> u64 {
            self.end
        }
    }

    impl<'a> fmt::Debug for Bucket<'a> {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "Bucket {{ {}..{} }}", self.start, self.end)
        }
    }
}
