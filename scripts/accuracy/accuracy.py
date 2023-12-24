# %%
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.utils import shuffle
import matplotlib.pyplot as plt
import seaborn as sns
from elyawy.io import load_sims_df
from elyawy.constants import SUMSTATS_LIST
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
# %%

_parser = argparse.ArgumentParser(allow_abbrev=False)

_parser.add_argument("-i", "--input", action='store',metavar="Input folder", type=str, required=True)
_parser.add_argument("-s", "--size", action='store',metavar="Threshold", type=int, required=True)
_parser.add_argument("-e", "--epsilon", action='store',metavar="Threshold", type=int, required=True)
_parser.add_argument('-a','--aligner', action='store',metavar="Alignment program to use" , type=str, required=True)

args = _parser.parse_args()

MAIN_PATH = Path(args.input)
SAMPLE_SIZE = args.size
EPSILON = args.epsilon
ALIGNER_NAME = args.aligner

res_path = MAIN_PATH.resolve()
print(res_path)
full_data, regressors, reg_stats = load_sims_df(data_path=res_path, correction=True, aligner=ALIGNER_NAME, sample_size=SAMPLE_SIZE)

# %%
full_zipf = shuffle(full_data[full_data["length_distribution"] == "zipf"])
full_geometric = shuffle(full_data[full_data["length_distribution"] == "geometric"])
full_poisson = shuffle(full_data[full_data["length_distribution"] == "poisson"])
num_test = 100
test_zipf = full_zipf[:num_test]
test_geometric = full_geometric[:num_test]
test_poisson = full_poisson[:num_test]

test_data = shuffle(pd.concat([test_zipf, test_geometric, test_poisson]))
remaining_data = shuffle(pd.concat([full_zipf[num_test:],full_geometric[num_test:], full_poisson[num_test:]]))

# %%
test_data_sum_stats = test_data[SUMSTATS_LIST].astype("float")
remaining_sum_stats = remaining_data[SUMSTATS_LIST].astype("float")
cov = np.cov(remaining_sum_stats.T)
cov = cov + np.eye(len(cov))*1e-4
inv_covmat = np.linalg.inv(cov)


# %%
predicted_dist = remaining_data["length_distribution"].reset_index(drop=True)
true_dists = test_data["length_distribution"].reset_index(drop=True)

# %%
abc_counts = []
count = 0
for u in test_data_sum_stats.values:
    def mahalanobis(simulated_stats):
        u_minus_v = simulated_stats-u
        left = np.dot(u_minus_v, inv_covmat)
        mahal = np.sqrt(np.sum(u_minus_v*left, axis=1))
        return mahal
    distances_from_u = mahalanobis(remaining_sum_stats).reset_index(drop=True)
    distances_from_u.name = "distances"
    distances_from_u = pd.concat([distances_from_u, predicted_dist], axis=1)
    distances_from_u.columns = ["distances", "predicted_distribution"]

    dist_distribution = distances_from_u.nsmallest(EPSILON, columns="distances")
    dist_distribution = dist_distribution["predicted_distribution"].value_counts()
    dist_distribution.name = f'{count}'
    count+=1
    abc_counts.append(dist_distribution)
all_counts = pd.concat(abc_counts,  axis=1)


# %%
abc_winners = pd.DataFrame(all_counts.fillna(0).idxmax())
abc_winners.columns = ["predicted"]

# %%
td_list = true_dists[:count]
td_list.name = "true_dists"
td_list = pd.DataFrame(td_list).T
td_list.columns = all_counts.columns
final_results = pd.concat([all_counts.fillna(0), abc_winners.T, td_list]).T

# %%
final_results.to_csv(res_path / f"simulation_results_{EPSILON}_{SAMPLE_SIZE}.csv")

# %%
confusion_matrix = pd.crosstab(final_results["true_dists"], final_results["predicted"])

# %%
sns.heatmap(confusion_matrix, annot=True, cmap='Blues', fmt='g')
plt.savefig(res_path / f"selection_results_{EPSILON}_{SAMPLE_SIZE}.png", dpi=200)


