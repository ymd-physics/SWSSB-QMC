use crate::{aux::{EMPTY_LEG, NULL_OP}, setup::SETUP, unsafe_get, unsafe_modify};

#[allow(non_snake_case)]
impl SETUP {
    fn randomly_freeze_Ws(&mut self) {
        for r in 0..self.Nr {
            unsafe {
                if unsafe_get!(self.rep_W_bs, r).is_some() {
                    unsafe_modify!(self.rep_W_frozen, r) = Some(self.rand_prob() <= 0.5);
                }
            }
        }
    }

    pub fn make_vertex_list(&mut self) {
        // ==================================================
        //  Initialization
        // ==================================================
        self.randomly_freeze_Ws();   // attention!
        
        for v in 0..4 * self.total_t {
            unsafe { unsafe_modify!(self.vertex_list, v) = EMPTY_LEG; }
        }

        for r in 0..self.Nr {
            for s in 0..self.N {
                unsafe { unsafe_modify!(self.rep_v_front, r, s) = EMPTY_LEG; }
                unsafe { unsafe_modify!(self.rep_v_back, r, s) = EMPTY_LEG; }
            }
        }

        // ==================================================
        //  Linking INSIDE the replicas        
        // ==================================================
        self.link_inside_rho0_reps();
        self.link_inside_mini_reps();

        // =========================================================================================
        //  Make the boundary condition
        //      Note that there is no free spin or EMPTY_LEG in in the 0th replica
        //      We use only the spins in rep-0 to track the cluster update
        // =========================================================================================
        for s in 0..self.N {
            unsafe {
                // sequentially connect the 0th replica with other replicas
                for r in 1..self.Nr {
                    let s_v_front = unsafe_get!(self.rep_v_front, r, s);
                    if s_v_front != EMPTY_LEG {
                        let s_v_back = unsafe_get!(self.rep_v_back, 0, s);
                        unsafe_modify!(self.vertex_list, s_v_back as usize) = s_v_front;
                        unsafe_modify!(self.vertex_list, s_v_front as usize) = s_v_back;
                        unsafe_modify!(self.rep_v_back, 0, s) = unsafe_get!(self.rep_v_back, r, s); 
                    }
                }

                // connect the front and back legs in the 0th replica
                let s_v_front = unsafe_get!(self.rep_v_front, 0, s);
                let s_v_back = unsafe_get!(self.rep_v_back, 0, s);
                unsafe_modify!(self.vertex_list, s_v_front as usize) = s_v_back;
                unsafe_modify!(self.vertex_list, s_v_back as usize) = s_v_front;
            }
        }
    }
}

#[allow(non_snake_case)]
impl SETUP {
    fn make_site_op_legs(&mut self, t: usize, r: usize, s: usize) {
        unsafe {
            let s_v_last = unsafe_get!(self.rep_v_back, r, s);
            let v_leg0 = (4 * t) as i32;

            if s_v_last != EMPTY_LEG {
                unsafe_modify!(self.vertex_list, s_v_last as usize) = v_leg0;
                unsafe_modify!(self.vertex_list, v_leg0 as usize) = s_v_last;
            } 

            else {
                unsafe_modify!(self.rep_v_front, r, s) = v_leg0;
            }
            unsafe_modify!(self.rep_v_back, r, s) = v_leg0 + 2;
        }   
    }

    fn make_bond_op_legs(&mut self, t: usize, r: usize, b:usize) {
        unsafe {
            let s0 = unsafe_get!(self.b_sites, b, 0);
            let s1 = unsafe_get!(self.b_sites, b, 1);
            let v_leg0 = (4 * t) as i32;
                    
            // ----------------------------------------------------
            let s0_v_last = unsafe_get!(self.rep_v_back, r, s0);
            if s0_v_last != EMPTY_LEG {
                unsafe_modify!(self.vertex_list, s0_v_last as usize) = v_leg0;
                unsafe_modify!(self.vertex_list, v_leg0 as usize) = s0_v_last;
            } 
            
            else {
                unsafe_modify!(self.rep_v_front, r, s0) = v_leg0;
            }
                
            unsafe_modify!(self.rep_v_back, r, s0) = v_leg0 + 2;
            
            // ----------------------------------------------------
            let s1_v_last = unsafe_get!(self.rep_v_back, r, s1);
            if s1_v_last != EMPTY_LEG {
                unsafe_modify!(self.vertex_list, s1_v_last as usize) = v_leg0 + 1;
                unsafe_modify!(self.vertex_list, (v_leg0 + 1) as usize) = s1_v_last;
            } 
            
            else {
                unsafe_modify!(self.rep_v_front, r, s1) = v_leg0 + 1;
            }
                
            unsafe_modify!(self.rep_v_back, r, s1) = v_leg0 + 3;
        }
    }

    fn link_inside_rho0_reps(&mut self) {
        for ell in 0..self.rho0_rep_labels.len() {
            unsafe {
                let r = unsafe_get!(self.rho0_rep_labels, ell);

                // ----------------------------------------------------------------------
                //  Link front site-Ids
                // ----------------------------------------------------------------------
                for t in unsafe_get!(self.rep_t_start, r)..(unsafe_get!(self.rep_t_start, r) + self.N) {
                    // assert!(self.op_string[t] == ID_SITE);
                    let s0 = t - self.rep_t_start[r];     // attention here!
                    // assert!(s0.within(0, self.N));
                    self.make_site_op_legs(t, r, s0);
                }

                // ----------------------------------------------------------------------
                //  Link Hamiltonian operators
                // ----------------------------------------------------------------------
                for t in (unsafe_get!(self.rep_t_start, r) + self.N)..(unsafe_get!(self.rep_t_start, r) + self.N + self.M) {
                    let op = unsafe_get!(self.op_string, t);
                    // assert!(op != ID_BOND && op != ID_SITE && op != XX_BOND);

                    if op != NULL_OP {
                        //  Bond operator
                        if op % 4 == 2 {
                            let b_t = (op / 4) as usize;
                            self.make_bond_op_legs(t, r, b_t);
                        }

                        //  Site operator
                        else {
                            let s0 = (op / 4) as usize;
                            self.make_site_op_legs(t, r, s0);
                        }
                    }
                }

                // ----------------------------------------------------------------------
                //  Link back site-Ids
                // ----------------------------------------------------------------------
                for t in (unsafe_get!(self.rep_t_end, r) - self.N)..unsafe_get!(self.rep_t_end, r) {
                    // assert!(self.rep_t_end[r] - self.N == self.rep_t_start[r] + self.N + self.M);
                    // assert!(self.op_string[t] == ID_SITE);

                    let s0 = t - (unsafe_get!(self.rep_total_t, r) - self.N) - unsafe_get!(self.rep_t_start, r);     // attention here!
                    // assert!(s0.within(0, self.N));
                    self.make_site_op_legs(t, r, s0);
                }
            }
        }
    }

    fn link_inside_mini_reps(&mut self) {
        for ell in 0..self.mini_rep_labels.len() {
            unsafe {
                let r = unsafe_get!(self.mini_rep_labels, ell);

                if let Some(b_W) = unsafe_get!(self.rep_W_bs, r) {
                    let s0 = unsafe_get!(self.b_sites, b_W, 0);
                    let s1 = unsafe_get!(self.b_sites, b_W, 1);

                    // ================================================================================
                    //  Link in this mini-rep with W-operators
                    // ================================================================================
                    if unsafe_get!(self.rep_W_frozen, r).unwrap() == false {
                        let t = unsafe_get!(self.rep_t_start, r);     // attention here! the zero time!
                        // assert!(self.op_string[t] == XX_BOND || self.op_string[t] == ID_BOND);
                        self.make_bond_op_legs(t, r, b_W);
                    }

                    self.make_site_op_legs(unsafe_get!(self.rep_t_start, r) + 1, r, s0);
                    self.make_site_op_legs(unsafe_get!(self.rep_t_start, r) + 2, r, s1);
                    self.make_site_op_legs(unsafe_get!(self.rep_t_start, r) + 3, r, s0);
                    self.make_site_op_legs(unsafe_get!(self.rep_t_start, r) + 4, r, s1);

                    if unsafe_get!(self.rep_W_frozen, r).unwrap() == false {
                        let t = unsafe_get!(self.rep_t_end, r) - 1;     // attention here! the last time
                        // assert!(self.rep_t_end[r] - 1 == self.rep_t_start[r] + 5);
                        // assert!(self.op_string[t] == XX_BOND || self.op_string[t] == ID_BOND);
                        self.make_bond_op_legs(t, r, b_W);
                    }

                    // ================================================================================
                    //  Link in the corresponding mini-rep WITHOUT W-operators
                    // ================================================================================
                    let corrspd_r = self.Nr - r;

                    self.make_site_op_legs(unsafe_get!(self.rep_t_start, corrspd_r) + 1, corrspd_r, s0);
                    self.make_site_op_legs(unsafe_get!(self.rep_t_start, corrspd_r) + 2, corrspd_r, s1);
                    self.make_site_op_legs(unsafe_get!(self.rep_t_start, corrspd_r) + 3, corrspd_r, s0);
                    self.make_site_op_legs(unsafe_get!(self.rep_t_start, corrspd_r) + 4, corrspd_r, s1);
                }
            }
        }
    }
}