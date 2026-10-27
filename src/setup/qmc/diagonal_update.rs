/* ==================================================================
    I consider op = 4 * b + a
    then
        op := -1, null operator
        op := 4 * s + 0, diag site operator
        op := 4 * s + 1, off-diag site operator
        op := 4 * b + 2, diag bond operator 
================================================================== */
use crate::{aux::{ID_BOND, ID_SITE, NULL_OP, XX_BOND}, setup::SETUP, unsafe_get, unsafe_modify};

impl SETUP {
    pub fn diagonal_update(&mut self) {
        for r in 0..self.Nr {
            // ----------------------------------------------------------------------
            //  We only need to consider the two copies of rho0
            // ----------------------------------------------------------------------
            if self.rho0_rep_labels.contains(&r) {
                unsafe {
                    for t in (unsafe_get!(self.rep_t_start, r) + self.N)..(unsafe_get!(self.rep_t_end, r) - self.N) {
                        // skip the site-identities
                        let op = unsafe_get!(self.op_string, t);

                        match op {
                            ID_SITE | ID_BOND | XX_BOND => { unreachable!("No site identity or W operator is expected."); }
                            
                            NULL_OP => {
                                if self.rand_prob() * ((self.M - unsafe_get!(self.rep_ns, r)) as f64) <= self.add_factor {
                                    if self.rand_prob() < self.selection_prob {
                                        let new_site = self.rand_site();
                                        unsafe_modify!(self.op_string, t) = (4 * new_site) as i32;
                                        unsafe_modify!(self.rep_ns, r) += 1;
                                    }

                                    else {
                                        let new_bond = self.rand_bond();
                                        let sp0 = unsafe_get!(self.rep_spins, 0, unsafe_get!(self.b_sites, new_bond, 0));
                                        let sp1 = unsafe_get!(self.rep_spins, 0, unsafe_get!(self.b_sites, new_bond, 1));
                                        if sp0 == sp1 {
                                            unsafe_modify!(self.op_string, t) = (4 * new_bond + 2) as i32;
                                            unsafe_modify!(self.rep_ns, r) += 1;
                                        }
                                    }
                                }
                            }

                            _ => {
                                let op_type = op % 4;
                                if op_type == 0 {
                                    if self.rand_prob() <= self.remove_factor * ((self.M - unsafe_get!(self.rep_ns, r) + 1) as f64) {
                                        unsafe_modify!(self.op_string, t) = NULL_OP;
                                        unsafe_modify!(self.rep_ns, r) -= 1;
                                    }
                                }

                                else if op_type == 2 {
                                    if self.rand_prob() <= self.remove_factor * ((self.M - unsafe_get!(self.rep_ns, r) + 1) as f64) {
                                        unsafe_modify!(self.op_string, t) = NULL_OP;
                                        unsafe_modify!(self.rep_ns, r) -= 1;
                                    }
                                }

                                else {
                                    let the_site = (op / 4) as usize;
                                    unsafe_modify!(self.rep_spins, 0, the_site) *= -1;
                                }
                            }
                        }
                    }
                }
            }

            // ----------------------------------------
            //  Update the rep-spins
            // ----------------------------------------
            if r != self.Nr - 1 {
                for s in 0..self.N {
                    unsafe { unsafe_modify!(self.rep_spins, r + 1, s) = unsafe_get!(self.rep_spins, 0, s); }
                }
            }
        }
    }
}