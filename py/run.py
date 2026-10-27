import numpy as np 
import os 
import time 

nohup = True

# ====================================
#   Setting params
# ====================================
Lx = 12
Ly = 1
beta = 8
J = 0.123
h = 2.345

#p_list_f64 = np.arange(start=0.35, stop=0.4 + 1e-9, step=0.005)
# p_list_f64 = np.arange(start=0.05, stop=0.14 + 1e-9, step=0.01)
p_list_f64 = np.array([0.3])

p_list = []
for p in p_list_f64:
    p_list.append("%.5f" % p)
p_list = np.array(p_list)

N_thm = 20000
N_rethm = 1000
N_stat = 500000

N_cores = 14
N_bins = 10

# ====================================
#   Generating PRNG seeds
# ====================================
seed0 = 77                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
# seed0 = int(time.time())
seed_list = [seed0]
N_seeds = len(p_list) * N_cores
for ell in range(1, N_seeds):
    seed_list.append(seed0 + 7 * ell)
seed_list = np.array(seed_list)

# ====================================
#   Initialization
# ====================================
os.system("cargo build --release")
os.system("rm -rf data && mkdir data")
os.system("rm -rf logs && mkdir logs")
time.sleep(1)   # to ensure that the target exe has been successfully built

path = f"./data/Lx{Lx}_Ly{Ly}_beta{beta}_J{J}_h{h}"
os.system(f"mkdir {path}")
np.savetxt("./path.txt", [path], fmt="%s", newline="")

# ==========================================
#   Save other params
# ==========================================
os.system(f"mkdir {path}/params")
np.savetxt(f"{path}/params/p_list.txt", p_list, fmt="%s")
np.savetxt(f"{path}/params/seed_list.txt", seed_list, fmt="%d")
np.savetxt(f"{path}/params/N_cores.txt", [N_cores], fmt="%d", newline="")
np.savetxt(f"{path}/params/N_bins.txt", [N_bins], fmt="%d", newline="")
np.savetxt(f"{path}/params/N_thm.txt", [N_thm], fmt="%d", newline="")
np.savetxt(f"{path}/params/N_rethm.txt", [N_rethm], fmt="%d", newline="")
np.savetxt(f"{path}/params/N_stat.txt", [N_stat], fmt="%d", newline="")

# ==========================================
#   Run the program
# ==========================================
print("-" * 100)
ell = 0

for p in p_list:
    print(f"♦︎ Submitting for p = {p}...")
    os.system(f"mkdir {path}/p{p}")

    for k in range(N_cores):
        print(f"\t--> Job {k}...")
        seed = seed_list[ell]   # get the seed
        subpath = f"{path}/p{p}/core{k}"
        os.system(f"mkdir {subpath}")

        all_params = f"{Lx} {Ly} {beta} {J} {h} {p} {N_thm} {N_rethm} {N_stat} {N_bins} {seed} {subpath}"
        command = f"./target/release/swssb_1d_unsafe {all_params}"

        if nohup:
            os.system(f"nohup {command} > ./logs/p{p}_core{k}.out 2>&1 &")
        else:
            os.system(f"{command}")

        time.sleep(0.2)
        ell += 1

print("♦︎ Done.")

# ==========================================
#  Print info
# ==========================================
print("-" * 100)
print(f"♦︎ Path = {path}")
print(f"♦︎ N_thm = {N_thm}, N_rethm = {N_rethm}, N_stat = {N_stat}")
print(f"♦︎ N_bins in each core = {N_bins}, N_cores for each p = {N_cores}")
print(f"♦︎ Total number of cores = {len(p_list) * N_cores}")
print(f"♦︎ Number of p points = {len(p_list)}")
print("-" * 100)
