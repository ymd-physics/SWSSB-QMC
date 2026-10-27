/****************************************************************
*  QMC simulations on Rényi-2 correlator
*  Author: Yi-Ming Ding
*  Email: dingyiming@westlake.edu.cn
*  Created: Dec 15, 2025
*  Updated: Dec 20, 2025
***************************************************************/
mod aux;
mod setup;
use std::{env, time::Instant, io::Write};

#[allow(non_snake_case)]
fn main() {
    let start_time = Instant::now();

    // ================================
    //  Collect params from shell
    // ================================
    let args: Vec<String> = env::args().collect();
    let param_Lx: usize = args[1].parse().unwrap();
    let param_Ly: usize = args[2].parse().unwrap();
    let param_beta: f64 = args[3].parse().unwrap();
    let param_J: f64 = args[4].parse().unwrap();
    let param_h: f64 = args[5].parse().unwrap();
    let param_p: f64 = args[6].parse().unwrap();

    let N_thm: usize = args[7].parse().unwrap();
    let N_rethm: usize = args[8].parse().unwrap();
    let N_stat: usize = args[9].parse().unwrap();
    let N_bins: usize = args[10].parse().unwrap();

    let param_seed: u64 = args[11].parse().unwrap();
    let path: String = args[12].parse().unwrap();

    // ================================================================
    //  Report the environment
    // ================================================================
    println!("■ QMC simulations on Rényi-2 correlator, with rho0 prepared from TFIM");
    println!("■ {param_Lx} x {param_Ly} lattices, beta = {param_beta}, J = {param_J}, h = {param_h}, p = {param_p}, seed = {param_seed}");
    println!("■ N_thm = {N_thm}, N_rethm = {N_rethm}, N_stat = {N_stat}, N_bins = {N_bins}");  

    // ====================================
    //  Preparing for writing the results
    // ====================================
    let mut file_renyi2_zz = aux::create_new_file(format!("{}/renyi2_zz.txt", path));
    let mut file_renyi2_mag = aux::create_new_file(format!("{}/renyi2_mag.txt", path));
    let mut file_renyi2_mag2 = aux::create_new_file(format!("{}/renyi2_mag2.txt", path));
    let mut file_renyi2_mag4 = aux::create_new_file(format!("{}/renyi2_mag4.txt", path));
    let mut file_renyi2_br = aux::create_new_file(format!("{}/renyi2_br.txt", path));

    // ===========================
    //  Monte Carlo simulations
    // ===========================
    let mut model = setup::SETUP::new(param_Lx, param_Ly, param_beta, param_J, param_h, param_seed, param_p);
    model.init();

    println!("\t---> Making rho0...");
    for _ in 0..N_thm {
        model.make_rho0();
    }
    println!("\t\t* MAX cut-off for rho0: {}", model.M);
    
    println!("\t---> Initializing for rho1...");
    model.init_reps();
    println!("\t\t* Total number of replicas: {}", model.Nr);
    println!("\t\t* Total length of rho1 x rho1: {}", model.total_t);
    
    println!("\t---> Rethermalizing...");
    for _ in 0..N_rethm {
        model.mc_update();
    }

    println!("\t---> Measuring...");
    for b in 0..N_bins {
        println!("\t\t# bin {}", b);

        model.init_measure_1d();
        for _ in 0..N_stat {
            model.mc_update();
            model.measure_1d();
        }
        model.statisticize_1d(N_stat as f64);

        // save the data 
        file_renyi2_zz.write_all(format!("{}\n", model.renyi2_zz_corr.iter().map(|x| format!("{:<20.15}", x)).collect::<Vec<_>>().join("\t")).as_bytes()).unwrap();
        file_renyi2_mag.write_all(format!("{:<20.15}\n", model.renyi2_mag).as_bytes()).unwrap();
        file_renyi2_mag2.write_all(format!("{:<20.15}\n", model.renyi2_mag2).as_bytes()).unwrap();
        file_renyi2_mag4.write_all(format!("{:<20.15}\n", model.renyi2_mag4).as_bytes()).unwrap();
        file_renyi2_br.write_all(format!("{:<20.15}\n", model.renyi2_br).as_bytes()).unwrap();
    }

    // ===========================
    //  Report the runtime
    // ===========================
    for _ in 0..50 { print!("-"); } print!("\n");
    aux::report_time(start_time);
}