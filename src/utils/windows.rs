// based on the slice::windows module in std
pub struct SeqWindows<'a> {
    pub seq: &'a [u8],
    pub size: usize,
    pub step: usize,
}

impl<'a> SeqWindows<'a> {
    pub fn new(slice: &'a [u8], size: usize, step: usize) -> Self {
        Self {
            seq: slice,
            size,
            step,
        }
    }
}

impl<'a> Iterator for SeqWindows<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<&'a [u8]> {
        if self.size > self.seq.len() {
            None
        } else {
            // get a slice from zero to window size
            let ret = Some(&self.seq[..self.size]);
            // modify the sequence so it starts the next iteration at
            // the step
            self.seq = &self.seq[self.step..];
            ret
        }
    }
}
