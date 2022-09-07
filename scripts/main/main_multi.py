import os, argparse, pathlib
from multiprocessing import Pool
from itertools import product
import numpy as np
import pandas as pd

from elyawy.constants import SUMSTATS_LIST, PARAMS_LIST, indel_models, length_distributions
from elyawy.sparta import Simulator, Msa

from sim_creator import SimConfig

ALL_MODELS = list(product(length_distributions, indel_models))

class P_Args:
    def __init__(self, p_input, p_numsim, p_lengthdist, p_model):
        self.input = p_input
        self.numsim = p_numsim
        self.lengthdist = p_lengthdist
        self.model = p_model



def generate_all_sims(args):
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
    for params in sim_params:
        numeric_params = [params[0],params[1], params[2], params[4], params[5]]

        sim.init_sim(*params)
        sim_msa = sim()

        sum_stats_all.append(np.array(numeric_params + sim_msa.get_sum_stats()))

    sum_stats_all = np.array(sum_stats_all)

    # check if correction exists and apply it

    data_full = np.concatenate([sum_stats_all, sim_params[:,[3]]], axis=1)
    data_full = pd.DataFrame(data_full, columns=PARAMS_LIST + SUMSTATS_LIST + ["length_distribution"])

    data_full.to_pickle(path_joiner(f"full_data_{LENGTH_DIST}_{INDEL_MODEL}.pkl"), compression="bz2")

    return "done"



if __name__ == "__main__":
    gr_parser = argparse.ArgumentParser(allow_abbrev=False)
    gr_parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)

    main_args = gr_parser.parse_args()

    DATA_PATH = pathlib.Path(main_args.input).resolve()

    print(DATA_PATH)

    processes_list = []
    for dir_path in DATA_PATH.iterdir():
        for model in ALL_MODELS:
            p_args = P_Args(str(dir_path), 100000, model[0], model[1])
            processes_list.append(p_args)


    with Pool(10) as p:
        print(p.map(generate_all_sims, processes_list))
