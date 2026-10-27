use crate::{aux::{EMPTY_LEG, FLIPPED, NULL_OP, UNFLIPPED}, setup::SETUP, unsafe_get, unsafe_modify};

#[allow(non_snake_case)]
impl SETUP {
    pub fn make_rho0(&mut self) {
        self.std_diagonal_update();
        self.std_cluster_update();
        self.adjust_m();
    }

    fn adjust_m(&mut self) {
        let new_M = self.n + self.n / 3;
        //let new_M = 2 * self.n;

        if self.M < new_M {
            self.op_string.extend(vec![NULL_OP; new_M - self.M]);
            self.M = new_M;
            self.vertex_list = vec![EMPTY_LEG; 4 * new_M];
            self.stack = Vec::with_capacity(8 * new_M);
        }
    }
}

/* ==================================================================
    I consider op = 4 * b + a
    then
        op := -1, null operator
        op := 4 * s + 0, diag site operator
        op := 4 * s + 1, off-diag site operator
        op := 4 * b + 2, diag bond operator 
================================================================== */
impl SETUP {
    fn std_diagonal_update(&mut self) {
        unsafe {
            for t in 0..self.M {
                let op = unsafe_get!(self.op_string, t);

                if op == NULL_OP {
                    // First decide whether to insert, then decide what to insert
                    if self.rand_prob() * ((self.M - self.n) as f64) <= self.add_factor {
                        if self.rand_prob() < self.selection_prob {
                            let new_site = self.rand_site();
                            unsafe_modify!(self.op_string, t) = (4 * new_site) as i32;
                            self.n += 1;
                        }

                        else {
                            let new_bond = self.rand_bond();
                            if unsafe_get!(self.spins, unsafe_get!(self.b_sites, new_bond, 0)) == unsafe_get!(self.spins, unsafe_get!(self.b_sites, new_bond, 1)) {
                                unsafe_modify!(self.op_string, t) = (4 * new_bond + 2) as i32;
                                self.n += 1;
                            }
                        }
                    }
                }

                else {
                    let op_type = op % 4;
                    if op_type == 0 || op_type == 2 {
                        if self.rand_prob() <= self.remove_factor * ((self.M - self.n + 1) as f64) {
                            unsafe_modify!(self.op_string, t) = NULL_OP;
                            self.n -= 1;
                        }
                    }

                    else {
                        let the_site = (op / 4) as usize;
                        unsafe_modify!(self.spins, the_site) *= -1;
                    }
                }
            }
        }
    }

    fn std_make_vertex_list(&mut self) {
        //  Initialize "vertex_list", "v_front", "v_back"
        for v in 0..4 * self.M {
            unsafe { unsafe_modify!(self.vertex_list, v) = EMPTY_LEG; }
        }
            
        for s in 0..self.N {
            unsafe {
                unsafe_modify!(self.v_front, s) = EMPTY_LEG;
                unsafe_modify!(self.v_back, s) = EMPTY_LEG;
            }
        }

        //  Make the vertex_list
        for t in 0..self.M {
            let op;
            unsafe { op = unsafe_get!(self.op_string, t); } 

            if op != NULL_OP {
                let v_leg0 = (4 * t) as i32;

                unsafe {
                    // ------------------------------
                    //  Bond operator
                    // ------------------------------
                    if op % 4 == 2 {
                        let b_t = (op / 4) as usize;
                        let s0 = unsafe_get!(self.b_sites, b_t, 0);
                        let s1 = unsafe_get!(self.b_sites, b_t, 1);
                        
                        {   
                            let s0_v_last = unsafe_get!(self.v_back, s0);
                            if s0_v_last != EMPTY_LEG {
                                unsafe_modify!(self.vertex_list, s0_v_last as usize) = v_leg0;
                                unsafe_modify!(self.vertex_list, v_leg0 as usize) = s0_v_last;
                            } 
                            
                            else {
                                unsafe_modify!(self.v_front, s0) = v_leg0;
                            }
                            
                            unsafe_modify!(self.v_back, s0) = v_leg0 + 2;
                        }
                        
                        {   
                            let s1_v_last = unsafe_get!(self.v_back, s1);
                            if s1_v_last != EMPTY_LEG {
                                unsafe_modify!(self.vertex_list, s1_v_last as usize) = v_leg0 + 1;
                                unsafe_modify!(self.vertex_list, (v_leg0 + 1) as usize) = s1_v_last;
                            } 
                            
                            else {
                                unsafe_modify!(self.v_front, s1) = v_leg0 + 1;
                            }
                            
                            unsafe_modify!(self.v_back, s1) = v_leg0 + 3;
                        }
                    }

                    // ------------------------------
                    //  Site operator
                    // ------------------------------
                    else {
                        let s0 = (op / 4) as usize;
                        let s0_v_last = unsafe_get!(self.v_back, s0);

                        if s0_v_last != EMPTY_LEG {
                            unsafe_modify!(self.vertex_list, s0_v_last as usize) = v_leg0;
                            unsafe_modify!(self.vertex_list, v_leg0 as usize) = s0_v_last;
                        } 

                        else {
                            unsafe_modify!(self.v_front, s0) = v_leg0;
                        }
                        
                        unsafe_modify!(self.v_back, s0) = v_leg0 + 2;
                    }
                }
            }
        }

        //  Make PBC correction
        let mut s_v_first;
        let mut s_v_last;

        unsafe {
            for s in 0..self.N {
                s_v_first = unsafe_get!(self.v_front, s);
                
                if s_v_first != -1 {
                    s_v_last = unsafe_get!(self.v_back, s);
                    unsafe_modify!(self.vertex_list, s_v_first as usize) = s_v_last;
                    unsafe_modify!(self.vertex_list, s_v_last as usize) = s_v_first;
                }
            }
        }
    }

    fn std_cluster_update(&mut self) {
        self.std_make_vertex_list();
        self.stack.clear();

        //  Make and flip the cluster
        for v in (0..4 * self.M).step_by(2) {
            unsafe {
                if unsafe_get!(self.vertex_list, v) < 0 {
                    // processed or with no operator
                    continue;
                }
            }

            self.flip = if self.rand_prob() <= 0.5 { FLIPPED } else { UNFLIPPED };
            self.stack.push(v);
            loop {
                if self.stack.is_empty() {
                    break;
                }
                self.std_make_cluster();
            }
        }

        //  Update the spins
        unsafe {
            for s in 0..self.N {
                if unsafe_get!(self.v_front, s) != EMPTY_LEG {
                    if unsafe_get!(self.vertex_list, unsafe_get!(self.v_front, s) as usize) == FLIPPED {
                        unsafe_modify!(self.spins, s) *= -1;
                    }
                }

                else if self.rand_prob() <= 0.5 {
                    unsafe_modify!(self.spins, s) *= -1;
                }
            }
        }
    }

    fn std_make_cluster(&mut self) {
        let v_start = self.stack.pop().unwrap();
        let v1; 
        unsafe { v1 = unsafe_get!(self.vertex_list, v_start); }

        if v1 < 0 { return; }
        
        unsafe {
            if unsafe_get!(self.vertex_list, v1 as usize) != self.flip {
                self.stack.push(v1 as usize);
            }
        }

        let t = v_start / 4;
        let op;
        unsafe { op = unsafe_get!(self.op_string, t); }
        let op_type = op % 4;

        if op_type < 2 {
            if self.flip == FLIPPED {
                unsafe { unsafe_modify!(self.op_string, t) ^= 1; }
            }
        }
        
        else {
            let v2 = v_start ^ 1;
            let v3 = v_start ^ 2;
            let v4 = v3 ^ 1;

            unsafe { if unsafe_get!(self.vertex_list, v2) >= 0 { self.stack.push(v2); } }
            unsafe { if unsafe_get!(self.vertex_list, v3) >= 0 { self.stack.push(v3); } }
            unsafe { if unsafe_get!(self.vertex_list, v4) >= 0 { self.stack.push(v4); } }
        }

        //  Mark "v_start" as processed
        self.vertex_list[v_start] = self.flip;
    }
}