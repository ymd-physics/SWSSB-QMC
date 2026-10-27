use prng_mt::mt19937::MT19937_64;
pub mod init;
pub mod qmc;
pub mod measure;
pub mod random;

#[allow(non_snake_case)]
pub struct SETUP {
    // ====================================================
    //  Basics
    // ====================================================
    pub Lx: usize,
    pub Ly: usize,
    pub beta: f64,
    pub J: f64,         
    pub h: f64,       
    rng: MT19937_64,

    // ====================================================
    //  Lattice related 
    // ====================================================
    pub N: usize,             
    Nb: usize,             
    b_sites: Vec<Vec<usize>>,    

    // ====================================================
    //  Configuration related 
    // ====================================================
    pub M: usize,         
    n: usize,      
    flip: i32,  

    spins: Vec<i32>,
    op_string: Vec<i32>,
    v_front: Vec<i32>,
    v_back: Vec<i32>,
    vertex_list: Vec<i32>,
    stack: Vec<usize>,

    // ====================================================
    //  Frequently-used values
    // ====================================================
    selection_prob: f64,
    add_factor: f64,
    remove_factor: f64,
    
    // ====================================================
    //  Replicated configuration related
    // ====================================================
    p: f64,     // decoherence strength
    pub Nr: usize,
    pub total_t: usize,
    pub rep_total_t: Vec<usize>,
    rep_ns: Vec<usize>,

    rep_spins: Vec<Vec<i32>>,
    rep_t_start: Vec<usize>,
    rep_t_end: Vec<usize>,

    rep_v_front: Vec<Vec<i32>>,
    rep_v_back: Vec<Vec<i32>>,

    rep_W_frozen: Vec<Option<bool>>,

    front_W_branch_to_neighbor: Vec<Option<bool>>,
    back_W_branch_to_neighbor: Vec<Option<bool>>,

    // tmp_spins: Vec<i32>,

    rho0_rep_labels: Vec<usize>,
    mini_rep_labels: Vec<usize>,

    rep_W_bs: Vec<Option<usize>>,
    rep_Id_tunneling: Vec<Vec<Option<bool>>>,

    // ====================================================
    //  Frequently-used values (replica)
    // ====================================================
    half_Nr: usize,
    merging_factor: f64,
    splitting_factor: f64,

    // ====================================================
    //  Measurements related
    // ====================================================
    r_L: usize,     // the replica related to the left space
    r_R: usize,     // the replica related to the right space
    
    pub renyi2_zz_corr: Vec<f64>,
    pub renyi2_mag: f64,
    pub renyi2_mag2: f64,
    pub renyi2_mag4: f64,
    pub renyi2_br: f64,
}