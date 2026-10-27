use crate::setup::SETUP;
pub mod make_rho0;
pub mod init_reps;
pub mod make_vertices;
pub mod cluster_update;
pub mod diagonal_update;
pub mod merge_n_split;

impl SETUP {
    pub fn mc_update(&mut self) {
        self.diagonal_update();
        self.merging_n_splitting();
        self.cluster_update();
    }   
}