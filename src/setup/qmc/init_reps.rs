use crate::{aux::{EMPTY_LEG, ID_BOND, ID_SITE}, setup::SETUP};

#[allow(non_snake_case)]
impl SETUP {
    pub fn init_reps(&mut self) {
        // =======================================================================
        //  Basic info
        // ======================================================================= 
        self.Nr = 4 * self.Nb + 2; 
        self.half_Nr = self.Nr / 2;

        self.rho0_rep_labels = vec![0, self.half_Nr];   

        for r in 0..self.Nr {
            if !self.rho0_rep_labels.contains(&r) { self.mini_rep_labels.push(r); }
        }

        assert!(self.rho0_rep_labels.len() + self.mini_rep_labels.len() == self.Nr);

        // =======================================================
        //  Init the rep-spins
        // =======================================================  
        // self.tmp_spins = self.spins.clone();
        self.rep_spins = vec![self.spins.clone(); self.Nr];

        // =======================================================
        //  Record the number of operators in each replicas
        // =======================================================  
        // note that "M" still represents the cut-off of a standard copy of state 
        for r in 0..self.Nr {
            if self.rho0_rep_labels.contains(&r) {
                self.rep_total_t.push(self.M + 2 * self.N);
                self.rep_ns.push(self.n);
            }
            else {
                assert!(self.mini_rep_labels.contains(&r));
                self.rep_total_t.push(6);   // 2 + 2 * 2 = 6
                self.rep_ns.push(0);
            }
        }

        assert!(self.rep_total_t.len() == self.Nr);
        assert!(self.rep_ns.len() == self.Nr);
        self.total_t = self.rep_total_t.iter().sum();

        // =======================================================
        //  Init the full "op_string"
        // =======================================================
        let mut new_op_string: Vec<i32> = Vec::with_capacity(self.total_t);
        
        for r in 0..self.Nr {
            if self.rho0_rep_labels.contains(&r) {
                // a copy of rho0
                new_op_string.extend(vec![ID_SITE; self.N]);
                new_op_string.extend_from_slice(&self.op_string[..]);
                new_op_string.extend(vec![ID_SITE; self.N]);
                assert!(self.rep_total_t[r] == 2 * self.N + self.M);
            }

            else {
                assert!(self.mini_rep_labels.contains(&r));
                new_op_string.extend(vec![ID_BOND; 1]);
                new_op_string.extend(vec![ID_SITE; 2]);
                new_op_string.extend(vec![ID_SITE; 2]);
                new_op_string.extend(vec![ID_BOND; 1]);
                assert!(self.rep_total_t[r] == 6);
            }
        }

        self.op_string = new_op_string;
        assert!(self.op_string.len() == self.total_t);

        // =======================================================
        //  Record the start and end times of each replica
        // =======================================================
        self.rep_t_start = vec![0];
        for r in 1..self.Nr {
            self.rep_t_start.push(self.rep_total_t[..r].iter().sum());
        }

        for r in 0..(self.Nr - 1) {
            self.rep_t_end.push(self.rep_t_start[r + 1]);
        }
        self.rep_t_end.push(self.total_t);

        assert!(self.rep_t_start.len() == self.Nr);
        assert!(self.rep_t_end.len() == self.Nr);

        // =====================================================================
        //  Init for W operators
        // =====================================================================
        self.rep_W_bs = vec![None; self.Nr];
        for r in 1..(1 + self.Nb) {
            self.rep_W_bs[r] = Some(r - 1);
        }
        for r in (self.half_Nr + 1)..(self.half_Nr + self.Nb + 1) {
            self.rep_W_bs[r] = Some(r - self.half_Nr - 1);
        }
    
        self.rep_W_frozen = vec![None; self.Nr];
        self.front_W_branch_to_neighbor = vec![None; self.Nr];
        self.back_W_branch_to_neighbor = vec![None; self.Nr];

        // =====================================================================
        //  Init for tunneling aux-ids
        // =====================================================================
        self.rep_Id_tunneling = vec![vec![None; self.N]; self.Nr];

        for r in 0..self.Nr {
            // only "2 * Nb" of them are Some(...)
            if let Some(b_W) = self.rep_W_bs[r] {
                let s0 = self.b_sites[b_W][0];
                let s1 = self.b_sites[b_W][1];
                self.rep_Id_tunneling[r][s0] = Some(true);
                self.rep_Id_tunneling[r][s1] = Some(true);

                let corrspd_r = self.Nr - r;    // Nr = 4 * Nb + 2
                self.rep_Id_tunneling[corrspd_r][s0] = Some(true);
                self.rep_Id_tunneling[corrspd_r][s1] = Some(true);

            }
        }
        assert!(self.rep_Id_tunneling.len() == self.Nr);

        // =======================================================
        //  Init for other data structures
        // =======================================================
        self.vertex_list = vec![EMPTY_LEG; 4 * self.total_t];
        self.stack.reserve(8 * self.total_t);
        
        self.rep_v_front = vec![vec![EMPTY_LEG; self.N]; self.Nr];
        self.rep_v_back = vec![vec![EMPTY_LEG; self.N]; self.Nr];

        // =======================================================
        //  Destroy irrelevant data structures
        // =======================================================
        self.spins.clear();
        self.v_front.clear();
        self.v_back.clear();
    }
}