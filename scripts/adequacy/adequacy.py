import pathlib, argparse, os
import numpy as np
import pandas as pd
import Sparta
from elyawy.constants import SUMSTATS_LIST
from elyawy.io import load_sims_df
from elyawy.sparta import Msa, Simulator

_parser = argparse.ArgumentParser(allow_abbrev=False)
_parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
# _parser.add_argument('-c','--config', action='store',metavar="Simulation config" , type=str, required=True)
_parser.add_argument('-n','--numsim', action='store',metavar="Number of simulations" , type=int, required=True)
# _parser.add_argument('-s','--seed', action='store',metavar="Simulation config" , type=int, required=False)
_parser.add_argument('-l','--lengthdist', action='store',metavar="Simulation config" , type=str, required=True)


args = _parser.parse_args()

MAIN_PATH = pathlib.Path(args.input).resolve()
NUM_SIMS = args.numsim
LENGTH_DISTRIBUTION = args.lengthdist


TREE_PATH = None
MSA_PATH = None

if len( n := list(MAIN_PATH.glob("*.tree"))) == 1:
    TREE_PATH = n[0]

if len( n := list(MAIN_PATH.glob("*.fasta"))) == 1:
    MSA_PATH = n[0]

if TREE_PATH is None or MSA_PATH is None:
    print("no fasta or tree file")
    exit()

simulator = Simulator(TREE_PATH)
empirical_msa = Msa(str(MSA_PATH))
empirical_sum_stats = empirical_msa.get_sum_stats()

# cleanup and preprocess
full_data, regressors, regressors_stats = load_sims_df(MAIN_PATH, correction=True)
full_data_columns = list(full_data.columns)

data_types = {column:float  for column in full_data.columns}
data_types["length_distribution"] = str
data_types["indel_model"] = str
full_data = full_data.astype(data_types)

reorder_columns = [ 
                    "root_length",
                    "length_distribution",
                    "length_param_insertion",
                    "length_param_deletion",
                    "insertion_rate",
                    "deletion_rate",
                    "indel_model"
                  ]

# distances calculation:
true_msa_sum_stats = np.array(empirical_sum_stats)
kept_stats_indices = list(range(len(SUMSTATS_LIST)))

if regressors_stats is not None:
    kept_stats_indices = regressors_stats[regressors_stats['pearsonr'] > 0.85].index
    if len(kept_stats_indices) < 15:
        kept_stats_indices = sorted(regressors_stats.nlargest(15, 'pearsonr').index)
    SUMSTATS_LIST = [SUMSTATS_LIST[i] for i in kept_stats_indices]
    true_msa_sum_stats = true_msa_sum_stats[kept_stats_indices]
simulated_sum_stats = full_data[SUMSTATS_LIST].astype(np.float32)


cov = np.cov(simulated_sum_stats.T)
inv_covmat = np.linalg.inv(cov)
u_minus_v = true_msa_sum_stats-simulated_sum_stats
left = np.dot(u_minus_v, inv_covmat)
mahalanobis_distances = np.sqrt(np.sum(u_minus_v*left, axis=1))

full_data["distances"] = mahalanobis_distances


top_number = 50


all_sims_data = pd.DataFrame()

top_sims = full_data[full_data["length_distribution"] == LENGTH_DISTRIBUTION].nsmallest(top_number, "distances")

params_top_sims = top_sims.loc[:,reorder_columns]
params_top_sims["length_distribution"] = top_sims["length_distribution"]
params_top_sims["indel_model"] = top_sims["indel_model"]

params_top_sims = params_top_sims[reorder_columns]
dist_samples = params_top_sims.sample(NUM_SIMS, replace=True).reset_index(drop=True)

def simulate_samples(params):
    params_dict = {
        "root_length": int(params["root_length"]),
        "length_distribution": params["length_distribution"],
        "length_parameter_insertion": params["length_param_insertion"],
        "length_parameter_deletion": params["length_param_deletion"],
        "insertion_rate": params["insertion_rate"],
        "deletion_rate": params["deletion_rate"],
    }
    simulator.init_sim(**params_dict)
    sim_msa = simulator()
        
    temp_params = params.values.tolist()[:-1]    
    temp_params = [temp_params[0]] + temp_params[2:]
    temp_params += sim_msa.get_sum_stats()

    temp_params = np.array(temp_params).reshape(1,-1)
    chosen_reg = f"{LENGTH_DISTRIBUTION}_{params.values[-1]}"
    if regressors is not None:
        temp_params = np.array([regressor.predict(temp_params).T for regressor in regressors[chosen_reg]])
    else:
        temp_params = (temp_params).T

    temp_params = temp_params[kept_stats_indices]
    return pd.Series(temp_params.reshape(-1).tolist(), SUMSTATS_LIST)
sample_sum_stats = dist_samples.apply(simulate_samples, axis=1)
# sample_sum_stats["length_distribution"] = LENGTH_DISTRIBUTION
# all_sims_data = pd.concat([all_sims_data, sample_sum_stats], axis=0)


data_full = pd.concat([sample_sum_stats, dist_samples], axis=1)

data_full = data_full[full_data_columns[:5] + SUMSTATS_LIST + full_data_columns[-2:]]
print(data_full)


try:
    os.mkdir(MAIN_PATH / "adequacy")
except:
    print("adequacy folder exists already")

data_full.to_csv(MAIN_PATH / "adequacy" / f"{LENGTH_DISTRIBUTION}.csv")