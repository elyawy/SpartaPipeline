import pickle, pathlib
from itertools import product

import numpy as np
import pandas as pd

from constants import *

def load_sims_df(data_path=None, correction=False, sample_size=100000):
    simulations_df = []
    regressors_stats = {} if correction is True else None
    regressors_dict = {} if correction is True else None
    all_models = product(length_distributions, indel_models)
    for model in all_models:
        lendist,indel_model = model[0], model[1]
        if correction is True:
            regressors_path = pathlib.Path(data_path, "correction", f"regressors_{lendist}_{indel_model}")
            if not regressors_path.exists():
                print("regressors file not found")
                exit(1)
            with open(regressors_path, 'rb') as f:
                regressors = pickle.load(f)
                regressors_dict[f"{lendist}_{indel_model}"] = regressors
            regressors_score_path = pathlib.Path(data_path, "correction", f"regression_performance_{lendist}_{indel_model}.csv")
            regressors_stats = pd.read_csv(regressors_score_path)

        simulations_path = pathlib.Path(data_path, f"full_data_{lendist}_{indel_model}.pkl")
        if not simulations_path.exists():
            print("simulations data not found")
            exit(1)
        temp_df = pd.read_pickle(simulations_path, compression='bz2')
        temp_df["indel_model"] = indel_model
        if temp_df.shape[0] < sample_size:
            print("sample size larger than number of simulations available")
            exit(1)
        temp_df = temp_df.sample(n=sample_size).reset_index(drop=True)
        if correction is True:
            
            temp_data = temp_df[PARAMS_LIST + SUMSTATS_LIST].values
            temp_data = np.array([regressor.predict(temp_data).T for regressor in regressors])
            temp_data = pd.DataFrame(temp_data.T, columns=SUMSTATS_LIST)
            temp_df[SUMSTATS_LIST] = temp_data
        simulations_df.append(temp_df)
    simulations_df = pd.concat(simulations_df, axis=0).reset_index(drop=True)
    return simulations_df, regressors_dict, regressors_stats