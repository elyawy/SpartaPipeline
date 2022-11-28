import pathlib, argparse
import numpy as np
import pandas as pd
from tqdm import tqdm

from elyawy.constants import SUMSTATS_LIST, PARAMS_LIST
from elyawy.sparta import Simulator, Msa

from sim_creator import SimConfig

gr_parser = argparse.ArgumentParser(allow_abbrev=False)
gr_parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
gr_parser.add_argument('-n','--numsim', action='store',metavar="Number of simulations" , type=int, required=True)
gr_parser.add_argument('-l','--lengthdist', action='store',metavar="Length distribution" , type=str, required=True)
gr_parser.add_argument('-m','--model', action='store',metavar="Indel model" , type=str, required=True)
gr_parser.add_argument('-v','--verbose', action='store_true')
gr_parser.set_defaults(verbose=False)

args = gr_parser.parse_args()

VERBOSE = args.verbose
if VERBOSE:
    print("\u0332".join("verbose mode".upper()) + ":")


MAIN_PATH = pathlib.Path(args.input).resolve()

TREE_PATH = None
MSA_PATH = None
for file in MAIN_PATH.iterdir():
    if ".tree" in str(file):
        TREE_PATH = str(MAIN_PATH / file)
    if ".fasta" in str(file):
        MSA_PATH = str(MAIN_PATH / file)

if TREE_PATH is None or MSA_PATH is None:
    print("no fasta or tree file")
    exit()

if VERBOSE:
    print("\u0332".join("Tree path")+ ": " + TREE_PATH)
    print("\u0332".join("Msa path")+ ": " + MSA_PATH)
    print()



NUM_SIMS = args.numsim



true_msa = Msa(MSA_PATH)
sim = Simulator(TREE_PATH)

true_msa_sum_stats = np.array(true_msa.get_sum_stats())
seq_lengths_in_msa = [true_msa.get_shortest_seq_length(), true_msa.get_longest_seq_length()]


sim_config = SimConfig(len_dist=args.lengthdist,
                       seq_lengths=seq_lengths_in_msa,
                       indel_model=args.model)

SEED = int(sim_config.seed)
LENGTH_DIST = args.lengthdist
INDEL_MODEL = args.model
# np.random.seed(int(SEED))
# minimum number of simulations for correction is 32, 200+ recommended.
sim_params = sim_config.get_random_sim(NUM_SIMS)

sum_stats_all = []
for params in tqdm(sim_params) if VERBOSE else sim_params:
    numeric_params = [params[0],params[1], params[2], params[4], params[5]]

    sim.init_sim(*params)
    sim_msa = sim()

    sum_stats_all.append(np.array(numeric_params + sim_msa.get_sum_stats()))


sum_stats_all = np.array(sum_stats_all)

# check if correction exists and apply it

data_full = np.concatenate([sum_stats_all, sim_params[:,[3]]], axis=1)
data_full = pd.DataFrame(data_full, columns=PARAMS_LIST + SUMSTATS_LIST + ["length_distribution"])

data_full.to_pickle(MAIN_PATH / f"full_data_{LENGTH_DIST}_{INDEL_MODEL}.pkl", compression="bz2")
