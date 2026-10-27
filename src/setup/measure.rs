use crate::{setup::SETUP, unsafe_get, unsafe_modify};

#[allow(non_snake_case)]
impl SETUP {
    fn get_1d_renyi2_zz_corr(&self, i: usize, j: usize) -> f64 {
        // For 1D "PBC", make use of the translation symmetry
        let mut val = 0.0;
        for ell in 0..self.N {
            unsafe {
                val += (
                    unsafe_get!(self.rep_spins, self.r_L, (i + ell) % self.N) * 
                    unsafe_get!(self.rep_spins, self.r_L, (j + ell) % self.N) * 
                    unsafe_get!(self.rep_spins, self.r_R, (i + ell) % self.N) * 
                    unsafe_get!(self.rep_spins, self.r_R, (j + ell) % self.N)
                ) as f64;
            }
        }
        val 
    }
}

#[allow(non_snake_case)]
impl SETUP {
    pub fn init_measure_1d(&mut self) {
        for s in 0..self.N { 
            unsafe { unsafe_modify!(self.renyi2_zz_corr, s) = 0.0; }
        }
        self.renyi2_mag = 0.0;
        self.renyi2_mag2 = 0.0;
        self.renyi2_mag4 = 0.0;
    }

    pub fn measure_1d(&mut self) {
        let mut maglr = 0.0;
        
        for s in 0..self.N {
            unsafe {
                //  Rényi-2 correlator
                unsafe_modify!(self.renyi2_zz_corr, s) += self.get_1d_renyi2_zz_corr(0, s);
                // Rényi-2 magnetization
                maglr += (unsafe_get!(self.rep_spins, self.r_L, s) * unsafe_get!(self.rep_spins, self.r_R, s)) as f64;
            }
        }

        // Rényi-2 magnetization and its moments
        self.renyi2_mag += maglr;
        self.renyi2_mag2 += maglr.powf(2.0);
        self.renyi2_mag4 += maglr.powf(4.0);
    }

    pub fn statisticize_1d(&mut self, mc_steps: f64) {
        for s in 0..self.N {
            unsafe { unsafe_modify!(self.renyi2_zz_corr, s) /= mc_steps * (self.N as f64); }
        }

        let N_f64 = self.N as f64;

        self.renyi2_mag /= mc_steps * N_f64;                    // attention!
        self.renyi2_mag2 /= mc_steps * N_f64.powf(2.0);         // attention!
        self.renyi2_mag4 /= mc_steps * N_f64.powf(4.0);         // attention!
        self.renyi2_br = self.renyi2_mag4 / (self.renyi2_mag2 * self.renyi2_mag2);
    }
}


// use crate::{setup::SETUP, unsafe_get, unsafe_modify};

// #[allow(non_snake_case)]
// impl SETUP {
//     fn get_2d_axes_info(&self, s: usize) -> (usize, usize, usize, usize){
//         let x_pos = s % self.Lx;
//         let y_pos = s / self.Ly;
//         let x_first = self.Lx * y_pos;
//         let y_first = x_pos;
//         (x_pos, y_pos, x_first, y_first)
//     }

//     fn get_2d_renyi2_zz_corr(&self, i: usize, j: usize) -> f64 {
//         // For 2D "PBC", make use of the translation symmetry
//         let (i_x_pos, i_y_pos, i_x_first, i_y_first) = self.get_2d_axes_info(i);
//         let (j_x_pos, j_y_pos, j_x_first, j_y_first) = self.get_2d_axes_info(j);

//         let mut val = 0.0;

//         for dx in 0..self.Lx {
//             let i_shift = i_x_first + (i_x_pos + dx) % self.Lx;
//             let j_shift = j_x_first + (j_x_pos + dx) % self.Lx;

//             unsafe {
//                 val += (
//                     unsafe_get!(self.rep_spins, self.r_L, i_shift) * unsafe_get!(self.rep_spins, self.r_L, j_shift) * 
//                     unsafe_get!(self.rep_spins, self.r_R, i_shift) * unsafe_get!(self.rep_spins, self.r_R, j_shift)
//                 ) as f64;
//             }
//         }

//         for dy in 0..self.Ly {
//             let i_shift = i_y_first + ((i_y_pos + dy) % self.Ly) * self.Lx;
//             let j_shift = j_y_first + ((j_y_pos + dy) % self.Ly) * self.Lx;

//             unsafe {
//                 val += (
//                     unsafe_get!(self.rep_spins, self.r_L, i_shift) * unsafe_get!(self.rep_spins, self.r_L, j_shift) * 
//                     unsafe_get!(self.rep_spins, self.r_R, i_shift) * unsafe_get!(self.rep_spins, self.r_R, j_shift)
//                 ) as f64;
//             }   
//         }
//         val 
//     }
// }

// #[allow(non_snake_case)]
// impl SETUP {
//     pub fn init_measure_2d(&mut self) {
//         for ell in 0..self.Lx { 
//             unsafe {unsafe_modify!(self.renyi2_zz_corr, ell) = 0.0; }
//         }
//         self.renyi2_mag = 0.0;
//         self.renyi2_mag2 = 0.0;
//         self.renyi2_mag4 = 0.0;
//     }

//     pub fn measure_2d(&mut self) {
//         //  Rényi-2 correlator
//         for ell in 0..self.Lx {
//             let s = ell * self.Lx + ell;
//             unsafe { unsafe_modify!(self.renyi2_zz_corr, ell) += self.get_2d_renyi2_zz_corr(0, s); }
//         }

//         // Rényi-2 magnetization and its moments
//         let mut maglr = 0.0;
//         for s in 0..self.N {
//             unsafe { maglr += (unsafe_get!(self.rep_spins, self.r_L, s) * unsafe_get!(self.rep_spins, self.r_R, s)) as f64; }
//         }

//         self.renyi2_mag += maglr;
//         self.renyi2_mag2 += maglr.powf(2.0);
//         self.renyi2_mag4 += maglr.powf(4.0);
//     }

//     pub fn statisticize_2d(&mut self, mc_steps: f64) {
//         for ell in 0..self.Lx {
//             // attention!!! 
//             unsafe { unsafe_modify!(self.renyi2_zz_corr, ell) /= mc_steps * (self.Lx + self.Ly) as f64;  } 
//         }

//         let N_f64 = self.N as f64;

//         self.renyi2_mag /= mc_steps * N_f64;                    // attention!
//         self.renyi2_mag2 /= mc_steps * N_f64.powf(2.0);         // attention!
//         self.renyi2_mag4 /= mc_steps * N_f64.powf(4.0);         // attention!
//         self.renyi2_br = self.renyi2_mag4 / (self.renyi2_mag2 * self.renyi2_mag2);
//     }
// }