use prng_mt::MT19937_64;
use crate::setup::SETUP;
use crate::aux::{EMPTY_LEG, NULL_OP};

#[allow(non_snake_case)]
impl SETUP {
    pub fn new(
        param_Lx: usize, 
        param_Ly: usize, 
        param_beta: f64, 
        param_J: f64, 
        param_h: f64, 
        param_seed: u64,
        param_p: f64,
    ) -> Self {
        Self {
            // ====================================================
            //  Basics
            // ====================================================
            Lx: param_Lx,
            Ly: param_Ly,
            beta: param_beta,
            J: param_J,
            h: param_h,
            rng: MT19937_64::new(param_seed),

            // ====================================================
            //  Lattice related 
            // ====================================================
            N: param_Lx * param_Ly,
            Nb: 0,                 // init in "init()"
            b_sites: Vec::new(),          // init in "init()"

            // ====================================================
            //  Configuration related 
            // ====================================================
            M: 10,
            n: 0,
            flip: 0,

            spins: Vec::new(),           // init in "init()"
            op_string: Vec::new(),       // init in "init()"
            v_front: Vec::new(),         // init in "init()"
            v_back: Vec::new(),          // init in "init()"
            vertex_list: Vec::new(),     // init in "init()"
            stack: Vec::new(),           // init in "init()"

            // ====================================================
            //  Frequently-used values
            // ====================================================
            selection_prob: 0.0,    // init in "init()"
            add_factor: 0.0,        // init in "init()"
            remove_factor: 0.0,     // init in "init()"
            
            // ====================================================
            //  Replicated configuration related
            // ====================================================
            p: param_p,
            total_t: 0,             // init later
            Nr: 0,                  // init later
            rep_total_t: Vec::new(),     // init later
            rep_ns: Vec::new(),     // init later

            rep_spins: Vec::new(),    // init later
            rep_t_start: Vec::new(),    // init later 
            rep_t_end: Vec::new(),      // init later

            rep_v_front: Vec::new(),    // init later
            rep_v_back: Vec::new(),     // init later

            rep_W_frozen: Vec::new(),   // init later

            front_W_branch_to_neighbor: Vec::new(),   // init later
            back_W_branch_to_neighbor: Vec::new(),    // init later

            // tmp_spins: Vec::new(),              // init later
            rho0_rep_labels: Vec::new(),        // init later
            mini_rep_labels: Vec::new(),        // init later

            rep_W_bs: Vec::new(),           // init later
            rep_Id_tunneling: Vec::new(), // init later

            // ====================================================
            //  Frequently-used values (replica)
            // ====================================================
            half_Nr: 0,                 // init later
            merging_factor: 0.0,          // init in "init()"
            splitting_factor: 0.0,        // init in "init()"

            // ====================================================
            //  Measurements related
            // ====================================================
            r_L: 0,         // init in "init()"
            r_R: 0,         // init in "init()"

            renyi2_zz_corr: vec![0.0; param_Lx],
            renyi2_mag: 0.0,
            renyi2_mag2: 0.0,
            renyi2_mag4: 0.0,
            renyi2_br: 0.0,
        }
    }

    pub fn init(&mut self) {
        // ====================================================
        //  Make a lattice
        // ====================================================
        assert!(self.Lx != 1, "For 1D lattice, set Ly == 1.");

        //  One-dimensional case
        if self.Ly == 1 {
            self.Nb = self.N;
            for i in 0..self.Nb {
                self.b_sites.push(vec![i, (i + 1) % self.N]);
            }
        }
        //  Two-dimensional case
        else {
            assert!(self.Lx == self.Ly, 
                "We require Lx == Ly for measuring correlation functions in the diagonal axis of the lattice.");

            self.Nb = 2 * self.Lx * self.Ly;
            for j in 0..self.Ly {
                for i in 0..self.Lx {
                    self.b_sites.push(vec![
                        self.Lx * j + i,
                        self.Lx * j + (i + 1) % self.Lx,
                    ]);
                }
            }

            for j in 0..self.Ly {
                for i in 0..self.Lx {
                    self.b_sites.push(vec![
                        self.Lx * j + i,
                        (self.Lx * j + self.Lx + i) % (self.Lx * self.Ly),
                    ])
                }
            }
        }
        
        // ====================================================
        //  Frequently-used values
        // ====================================================
        self.selection_prob = self.N as f64 * self.h
            / (self.Nb as f64 * 2.0 * self.J + (self.N as f64) * self.h);
        self.add_factor = self.beta * (self.N as f64 * self.h + 2.0 * self.J * self.Nb as f64);
        self.remove_factor = 1.0 / self.add_factor;

        // ====================================================
        //  Configuration related 
        // ====================================================
        self.spins = (0..self.N).map(|_| if self.rand_prob() > 0.5 { 1 } else { -1 }).collect();
        self.op_string = vec![NULL_OP; self.M];
        self.v_front = vec![EMPTY_LEG; self.N];
        self.v_back = vec![EMPTY_LEG; self.N];
        self.vertex_list = vec![EMPTY_LEG; 4 * self.M];
        self.stack.reserve(8 * self.M);

        // ====================================================
        //  Frequently-used values (replicas)
        // ====================================================
        self.r_L = self.Nb + 1;
        self.r_R = 3 * self.Nb + 2;

        assert!(self.p <= 1.0 && self.p >= 0.0);
        match self.p {
            1.0 => {
                self.splitting_factor = 0.0;
                self.merging_factor = 1.0;
            }
            0.0 => {
                self.splitting_factor = 1.0;
                self.merging_factor = 0.0;
            }
            _ => {
                self.splitting_factor = (1.0 - self.p) / (2.0 * self.p);
                self.merging_factor = (2.0 * self.p) / (1.0 - self.p);
            }
        }
    }
}