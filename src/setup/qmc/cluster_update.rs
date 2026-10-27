use crate::{aux::{FLIPPED, ID_BOND, ID_SITE, RangeCheck, UNFLIPPED, XX_BOND}, setup::SETUP, to_back, to_neighbor, to_opposite, unsafe_get, unsafe_modify};

#[allow(non_snake_case)]
impl SETUP {
    fn randomize_aux_branches(&mut self) {
        for ell in 0..self.mini_rep_labels.len() {
            unsafe {
                let r = unsafe_get!(self.mini_rep_labels, ell);

                if let Some(false) = unsafe_get!(self.rep_W_frozen, r) {
                    unsafe_modify!(self.back_W_branch_to_neighbor, r) = Some(self.rand_prob() <= 0.5);
                    unsafe_modify!(self.front_W_branch_to_neighbor, r) = Some(self.rand_prob() <= 0.5);
                }

                else {
                    unsafe_modify!(self.back_W_branch_to_neighbor, r) = None;
                    unsafe_modify!(self.front_W_branch_to_neighbor, r) = None;
                }
            }
        }
    }

    fn get_rep_label(&mut self, t: usize) -> usize {
        for r in 0..self.Nr {
            unsafe { if t.within(unsafe_get!(self.rep_t_start, r), unsafe_get!(self.rep_t_end, r)) { return r; } }
        } unreachable!("t = {t} is out of the range of all time!");
    }

    pub fn cluster_update(&mut self) {
        self.make_vertex_list();
        self.randomize_aux_branches();  // after making the vertex list
        self.stack.clear();

        // ===============================================================
        //  Grow and flip the cluster
        // ===============================================================
        for v in (0..4 * self.total_t).step_by(2) {
            unsafe { if unsafe_get!(self.vertex_list, v) < 0 { continue; } }

            self.flip = if self.rand_prob() <= 0.5 { FLIPPED } else { UNFLIPPED };
            self.stack.push(v);
            loop {
                if self.stack.is_empty() { break; }
                self.make_cluster();
            }
        }

        // =======================================================
        //  Update the spins saved in the 0th replica
        // =======================================================
        for s in 0..self.N {
            unsafe {
                if unsafe_get!(self.vertex_list, unsafe_get!(self.rep_v_front, 0, s) as usize) == FLIPPED {
                    unsafe_modify!(self.rep_spins, 0, s) *= -1;
                }
            }
        }
    }
    
    fn make_cluster(&mut self) {
        let v0 = self.stack.pop().unwrap();
        let v1;
        unsafe { v1 = unsafe_get!(self.vertex_list, v0); }

        // return if it has been processed
        if v1 < 0 { return; }   

        // link to the next one 
        if self.vertex_list[v1 as usize] != self.flip {
            self.stack.push(v1 as usize);
        }

        let t = v0 / 4;
        let op;
        unsafe { op = unsafe_get!(self.op_string, t); }

        match op {
            // ===========================================
            //  IDENTITY SITE OPERATORS
            // ===========================================
            ID_SITE => {
                unsafe {
                    // -------------------------------------------
                    //  Put the back of it into the cluster
                    // -------------------------------------------
                    let v_back = to_back!(v0);
                    if unsafe_get!(self.vertex_list, v_back) >= 0 { 
                        self.stack.push(v_back); 
                    }

                    // ----------------------------------------------------
                    //  Try tunneling if it is in a mini replica
                    // ----------------------------------------------------
                    let r = self.get_rep_label(t);

                    if self.rho0_rep_labels.contains(&r) == false {
                        let corrspd_r = self.Nr - r;

                        let b_W = unsafe_get!(self.rep_W_bs, r).or(unsafe_get!(self.rep_W_bs, corrspd_r)).unwrap();
                        
                        // assert!(t.within(self.rep_t_start[r] + 1, self.rep_t_end[r] - 1));
                        
                        let this_t = t - unsafe_get!(self.rep_t_start, r);
                        let s = {
                            if this_t % 2 == 1 { unsafe_get!(self.b_sites, b_W, 0) }
                            else { unsafe_get!(self.b_sites, b_W, 1) }
                        };

                        // assert!(s.within(0, self.N));

                        if let Some(true) = unsafe_get!(self.rep_Id_tunneling, r, s) {
                            let target_t = unsafe_get!(self.rep_t_start, corrspd_r) + this_t;
                            let target_v = 4 * target_t;

                            // assert!(self.op_string[target_t] == ID_SITE);

                            if unsafe_get!(self.vertex_list, target_v) >= 0 {
                                self.stack.push(target_v);
                            }
                        }
                    }
                }
            }   

            // ==========================================
            //  W operators
            // ==========================================
            ID_BOND | XX_BOND => {
                unsafe {
                    let r = self.get_rep_label(t);

                    if let Some(false) = unsafe_get!(self.rep_W_frozen, r) {
                        // assert!(self.mini_rep_labels.contains(&r));
                        // assert!(self.rep_W_bs[r].is_some());
                        // assert!(t == self.rep_t_start[r] || t == self.rep_t_end[r] - 1);

                        // ::::::::::::::::::::::::::::::::::::::::
                        //  If it at time-zero
                        // ::::::::::::::::::::::::::::::::::::::::
                        if t == unsafe_get!(self.rep_t_start, r) {
                            let this_to_neighbor = unsafe_get!(self.front_W_branch_to_neighbor, r).unwrap();   

                            // branching 
                            let v_branch = {
                                if this_to_neighbor { to_neighbor!(v0) }
                                else { to_opposite!(v0) }
                            };

                            if unsafe_get!(self.vertex_list, v_branch) >= 0 { 
                                self.stack.push(v_branch);
                                // avoid reduplicative flip
                                if self.flip == FLIPPED { 
                                    unsafe_modify!(self.op_string, t) ^= 1; 
                                }
                            }
                            
                            // let target_t = self.rep_t_end[r] - 1;
                            let target_to_neighbor = unsafe_get!(self.back_W_branch_to_neighbor, r).unwrap();

                            match (this_to_neighbor, target_to_neighbor) {
                                (true, true) => {
                                    // in this case, we have to ensure the connectivity on the outer sides 
                                    if v0 % 4 < 2 {
                                        let dt = 5;
                                        // assert!(t + dt == target_t);

                                        let target_v = to_back!(v0 + 4 * dt);
                                        if unsafe_get!(self.vertex_list, target_v) >= 0 {
                                            self.stack.push(target_v);
                                        }
                                    }
                                }

                                (true, false) => {
                                    // in this case, we just need to get the back of this leg
                                    // since everything would be connected
                                    let target_v = to_back!(v0);
                                    if unsafe_get!(self.vertex_list, target_v) >= 0 {
                                        self.stack.push(target_v);
                                    }
                                }

                                _ => {}
                            }
                        }

                        // ::::::::::::::::::::::::::::::::::::::::
                        //  If it at time-last
                        // ::::::::::::::::::::::::::::::::::::::::
                        else {
                            // assert!(t == self.rep_t_end[r] - 1);

                            let this_to_neighbor = unsafe_get!(self.back_W_branch_to_neighbor, r).unwrap();   

                            // branching 
                            let v_branch = {
                                if this_to_neighbor { to_neighbor!(v0) }
                                else { to_opposite!(v0) }
                            };

                            if unsafe_get!(self.vertex_list, v_branch) >= 0 { 
                                self.stack.push(v_branch);
                                // avoid reduplicative flip
                                if self.flip == FLIPPED { 
                                    unsafe_modify!(self.op_string, t) ^= 1; 
                                }
                            }

                            // let target_t = self.rep_t_start[r];
                            let target_to_neighbor = unsafe_get!(self.front_W_branch_to_neighbor, r).unwrap();

                            match (target_to_neighbor, this_to_neighbor) {
                                (true, true) => {
                                    // in this case, we have to ensure the connectivity on the outer sides 
                                    if v0 % 4 >= 2 {
                                        let dt = 5;
                                        // assert!(t - dt == target_t);

                                        let target_v = to_back!(v0 - 4 * dt);
                                        if unsafe_get!(self.vertex_list, target_v) >= 0 {
                                            self.stack.push(target_v);
                                        }
                                    }
                                }

                                (false, true) => {
                                    // in this case, we just need to get the back of this leg
                                    // since everything would be connected
                                    let target_v = to_back!(v0);
                                    if unsafe_get!(self.vertex_list, target_v) >= 0 {
                                        self.stack.push(target_v);
                                    }
                                }

                                _ => {}
                            }
                        }

                    }
                }
            }

            // ==========================================
            //  HAMILTONIAN OPERATORS
            // ==========================================
            _ => {
                let op_type = op % 4;

                // flip this site operator
                if op_type < 2 {
                    if self.flip == FLIPPED {
                        unsafe { unsafe_modify!(self.op_string, t) ^= 1; }
                    }
                }

                else {
                    let v2 = v0 ^ 1;
                    let v3 = v0 ^ 2;
                    let v4 = v3 ^ 1;
                    unsafe { if unsafe_get!(self.vertex_list, v2) >= 0 { self.stack.push(v2); } }
                    unsafe { if unsafe_get!(self.vertex_list, v3) >= 0 { self.stack.push(v3); } }
                    unsafe { if unsafe_get!(self.vertex_list, v4) >= 0 { self.stack.push(v4); } }
                }
            }
        }

        //  Mark "v0" as processed
        unsafe { unsafe_modify!(self.vertex_list, v0) = self.flip; }
    }
}