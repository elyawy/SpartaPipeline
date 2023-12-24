import argparse, pathlib
import numpy as np
import pandas as pd
# my special imports
import elyawy.io as eio
from elyawy.constants import SUMSTATS_LIST
from elyawy.sparta import Msa
from sklearn.base import BaseEstimator, TransformerMixin# define the transformer

class StandardMemoryScaler(BaseEstimator, TransformerMixin):

    def __init__(self, epsilon=1e-4):
        self._epsilon = epsilon
        
    def fit(self, X, y = None):
        self._mean = X.mean()
        self._std = X.std()

        return self

    def transform(self, X):
        X = (X-self._mean)/(self._std+self._epsilon)
       
        return X


_parser = argparse.ArgumentParser(allow_abbrev=False)
_parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
_parser.add_argument('-s','--size', action='store',metavar="Sample size", type=int, required=False, default=100000)
_parser.add_argument('-nc','--no-correction', action='store_false')


args = _parser.parse_args()

MAIN_PATH = pathlib.Path(args.input)
MAIN_PATH = MAIN_PATH if MAIN_PATH.exists() else exit(1)

MSA_PATH = None
for file in MAIN_PATH.iterdir():
    if ".fasta" == file.suffix:
        MSA_PATH = file

if MSA_PATH is None:
    print("no fasta file")
    exit(1)

true_msa_sum_stats = Msa(str(MSA_PATH)).get_sum_stats()
true_msa_sum_stats = np.array(true_msa_sum_stats)

CORRECTED = args.no_correction
SAMPLE_SIZE = args.size
sims_df, regs_dict,stats_reg_score_df = eio.load_sims_df(MAIN_PATH, correction=CORRECTED, aligner=ALIGNER_NAME,sample_size=SAMPLE_SIZE)

CORRECTION_THRESHOLD = 0.0

if stats_reg_score_df is not None:
    kept_stats_indices = stats_reg_score_df[stats_reg_score_df['pearsonr'] > CORRECTION_THRESHOLD].index
    if CORRECTION_THRESHOLD <= 0.0:
        kept_stats_indices = stats_reg_score_df.index
    if len(kept_stats_indices) < 15:
        kept_stats_indices = stats_reg_score_df.nlargest(15, 'pearsonr').index
    SUMSTATS_LIST = [SUMSTATS_LIST[i] for i in kept_stats_indices]
    true_msa_sum_stats = true_msa_sum_stats[kept_stats_indices]
simulated_sum_stats = sims_df[SUMSTATS_LIST].astype(np.float32)

dist_method = "mahal"

if dist_method == "mahal":
    cov = np.cov(simulated_sum_stats.T)
    cov = cov + np.eye(len(cov))*1e-4
    inv_covmat = np.linalg.inv(cov)
    u_minus_v = true_msa_sum_stats-simulated_sum_stats
    left = np.dot(u_minus_v, inv_covmat)
    all_distances = np.sqrt(np.sum(u_minus_v*left, axis=1))
if dist_method == "euclid":
    weights = 1/(simulated_sum_stats.std(axis=0) + 0.001)
    all_distances = np.sum(weights*(simulated_sum_stats - true_msa_sum_stats)**2, axis=1)


sims_df["distances"] = all_distances

top_data = sims_df.nsmallest(1, "distances")

corrected_str = 'corrected' if CORRECTED else 'not_corrected'

file_name = f"top_params_{dist_method}_{SAMPLE_SIZE}_{corrected_str}.csv"
results_path = pathlib.Path(MAIN_PATH,file_name)
top_data.to_csv(results_path)



