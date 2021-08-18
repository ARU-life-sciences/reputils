use std::iter::Peekable;
pub struct SequentialCount<I>
where
    I: Iterator,
{
    pub iter: Peekable<I>,
}

impl<I> SequentialCount<I>
where
    I: Iterator,
{
    pub fn new(iter: I) -> Self {
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
