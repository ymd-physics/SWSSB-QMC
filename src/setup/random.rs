use crate::setup::SETUP;
const PRECISION: u64 = 1u64 << 53;

impl SETUP {
    pub fn rand_prob(&mut self) -> f64 {
        let val = self.rng.next() >> (64 - 53);
        val as f64 / PRECISION as f64
    }

    pub fn rand_bond(&mut self) -> usize {
        self.rng.next() as usize % self.Nb
    }

    pub fn rand_site(&mut self) -> usize { 
        self.rng.next() as usize % self.N 
    }
}