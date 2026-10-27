use crate::{aux::XX_BOND, setup::SETUP, unsafe_get, unsafe_modify};

impl SETUP {
    pub fn merging_n_splitting(&mut self) {
        for r in 0..self.Nr {
            unsafe {
                // only "2 * Nb" of them are Some(...)
                if let Some(b) = unsafe_get!(self.rep_W_bs, r) {
                    let s0 = unsafe_get!(self.b_sites, b, 0);
                    let s1 = unsafe_get!(self.b_sites, b, 1);

                    let corrspd_r = self.Nr - r;    // Nr = 4 * Nb + 2

                    let this_is_tunnelled = unsafe_get!(self.rep_Id_tunneling, r, s0).unwrap();
                    // assert!(this_is_tunnelled == self.rep_Id_tunneling[r][s1].unwrap());
                    // assert!(this_is_tunnelled == self.rep_Id_tunneling[corrspd_r][s0].unwrap());
                    // assert!(this_is_tunnelled == self.rep_Id_tunneling[corrspd_r][s1].unwrap());

                    if this_is_tunnelled {
                        // Check whether the corresponding states are consistent
                        if unsafe_get!(self.op_string, unsafe_get!(self.rep_t_start, r)) == XX_BOND {
                            assert!(unsafe_get!(self.rep_spins, r, s0) == -1 * unsafe_get!(self.rep_spins, corrspd_r, s0));   
                            assert!(unsafe_get!(self.rep_spins, r, s1) == -1 * unsafe_get!(self.rep_spins, corrspd_r, s1));   
                        }
                        else {
                            assert!(unsafe_get!(self.rep_spins, r, s0) == unsafe_get!(self.rep_spins, corrspd_r, s0));   
                            assert!(unsafe_get!(self.rep_spins, r, s1) == unsafe_get!(self.rep_spins, corrspd_r, s1));   
                        }

                        // Try splitting
                        if self.rand_prob() <= self.splitting_factor {
                            unsafe_modify!(self.rep_Id_tunneling, r, s0) = Some(false);
                            unsafe_modify!(self.rep_Id_tunneling, r, s1) = Some(false);
                            unsafe_modify!(self.rep_Id_tunneling, corrspd_r, s0) = Some(false);
                            unsafe_modify!(self.rep_Id_tunneling, corrspd_r, s1) = Some(false);
                        }
                    }

                    else {
                        let mut spins_are_consistent = false;

                        if unsafe_get!(self.op_string, unsafe_get!(self.rep_t_start, r)) == XX_BOND {
                            if (unsafe_get!(self.rep_spins, r, s0) == -1 * unsafe_get!(self.rep_spins, corrspd_r, s0)) && 
                                (unsafe_get!(self.rep_spins, r, s1) == -1 * unsafe_get!(self.rep_spins, corrspd_r, s1)) {
                                spins_are_consistent = true;
                            }
                        }
                        else {
                            if (unsafe_get!(self.rep_spins, r, s0) == unsafe_get!(self.rep_spins, corrspd_r, s0)) && 
                                (unsafe_get!(self.rep_spins, r, s1) == unsafe_get!(self.rep_spins, corrspd_r, s1)) {
                                spins_are_consistent = true;
                            }
                        }

                        if spins_are_consistent {
                            // Try merging
                            if self.rand_prob() <= self.merging_factor {
                                unsafe_modify!(self.rep_Id_tunneling, r, s0) = Some(true);
                                unsafe_modify!(self.rep_Id_tunneling, r, s1) = Some(true);
                                unsafe_modify!(self.rep_Id_tunneling, corrspd_r, s0) = Some(true);
                                unsafe_modify!(self.rep_Id_tunneling, corrspd_r, s1) = Some(true);
                            }
                        }
                    }
                }   
            }
        }
    }
}