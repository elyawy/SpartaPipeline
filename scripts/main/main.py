
import os, pickle, argparse
import numpy as np
import pandas as pd
from tqdm import tqdm

import Sparta
from sim_creator import SimConfig

gr_parser = argparse.ArgumentParser(allow_abbrev=False)
gr_parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
gr_parser.add_argument('-n','--numsim', action='store',metavar="Number of simulations" , type=int, required=True)
gr_parser.add_argument('-l','--lengthdist', action='store',metavar="Length distribution" , type=str, required=True)
gr_parser.add_argument('-m','--model', action='store',metavar="Indel model" , type=str, required=True)

args = gr_parser.parse_args()

MAIN_PATH = args.input
path_joiner = lambda x: os.path.join(MAIN_PATH, x)

TREE_PATH = None
MSA_PATH = None
for file in os.listdir(MAIN_PATH):
    if ".tree" in file:
        TREE_PATH = path_joiner(file)
    if ".fasta" in file:
        MSA_PATH = path_joiner(file)

if TREE_PATH is None or MSA_PATH is None:
    print("no fasta or tree file")
    exit()

NUM_SIMS = args.numsim
NUM_SIM_TOP = 100
PARAMS_LIST = ["root_length","length_param_insertion", "length_param_deletion", "insertion_rate", "deletion_rate"]
SUMSTATS_LIST = [f'SS_{i}' for i in range(0,27)]


true_msa = Sparta.Msa(MSA_PATH)
sim = Sparta.Sim(TREE_PATH)

true_msa.calc_stats()
true_msa_sum_stats = np.array(true_msa.get_sum_stats())
seq_lengths_in_msa = [true_msa.get_shortest_seq_length(), true_msa.get_longest_seq_length()]


sim_config = SimConfig(len_dist=args.lengthdist,
                       seq_lengths=seq_lengths_in_msa,
                       indel_model=args.model)

SEED = int(sim_config.seed)
LENGTH_DIST = sim_config.length_distribution
INDEL_MODEL = sim_config.indel_model
# np.random.seed(int(SEED))

# minimum number of simulations for correction is 32, 200+ recommended.
sim_params = sim_config.get_random_sim(NUM_SIMS)

sum_stats_all = []
for params in (sim_params):
    numeric_params = [params[0],params[2][0], params[3][0], params[4], params[5]]
    # sim.set_seed(int(SEED))
    # print(params)
    sim.init_sim(*params)
    sim_msa = sim.run_sim()
    # print("done wiht sim")
    sim_msa.calc_stats()
    sum_stats_all.append(np.array(numeric_params + sim_msa.get_sum_stats()))

sum_stats_all = np.array(sum_stats_all)

# check if correction exists and apply it

data_full = np.concatenate([sum_stats_all, sim_params[:,[1]]], axis=1)
data_full = pd.DataFrame(data_full, columns=PARAMS_LIST + SUMSTATS_LIST + ["length_distribution"])
data_full.to_pickle(path_joiner(f"full_data_{LENGTH_DIST}_{INDEL_MODEL}.pkl"), compression="bz2")
