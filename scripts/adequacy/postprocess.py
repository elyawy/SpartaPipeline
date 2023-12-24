import pathlib, copy, re, pickle
from itertools import product
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8')
import matplotlib.ticker as mtick
import seaborn as sns
from Bio import Phylo
from elyawy.constants import SUMSTATS_LIST, SUMSTATS_DEFINITION, length_distributions
from elyawy.sparta import Msa
import raxml_parser
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

REMOTE_PATH = pathlib.Path("/groups/pupko/elyawygoda/length_distributions/all_outputs/results").resolve() # sys.argv[2]
REMOTE_PATH.exists()

CLASSIFICATIONS = {}
CLASSIFICATIONS["eggnog"] = pd.read_csv(REMOTE_PATH / "eggnog_classifications.csv", index_col="index")["0"]
CLASSIFICATIONS["yeast"] = pd.read_csv(REMOTE_PATH / "yeast_classifications.csv", index_col="index")["0"]
# CLASSIFICATIONS["adequacy_sims"] = pd.read_csv(RESULTS_PATH / "adequacy_sims_classifications.csv", index_col="index")["0"]

def get_classification(current_path: pathlib.Path, corrected=True):
    return CLASSIFICATIONS[current_path.parent.stem][current_path.stem]


def get_percentile_summary(sims_df, empirical_df):
    to_df = {}
    for sum_stat in filter(lambda x: "SS" in x, sims_df.columns):
        percentile = stats.percentileofscore(sims_df[sum_stat], empirical_df[sum_stat], kind="weak")
        stat_def = SUMSTATS_DEFINITION[sum_stat].replace("\n", "_")#.replace("_"," ").title()
        to_df[stat_def] =  [sum_stat, empirical_df[sum_stat], percentile]
    
    percentile_summary =  pd.DataFrame(to_df,index=["stat_id","empirical","percentile"]).T
    percentile_summary["retained"] = (percentile_summary.percentile.between(2.5,97.5)).astype(int)

    return percentile_summary

def get_percentile_summary_pretty(sims_df, empirical_df):
    to_df = {}
    for sum_stat in filter(lambda x: "SS" in x, sims_df.columns):
        percentile = stats.percentileofscore(sims_df[sum_stat], empirical_df[sum_stat], kind="weak")
        stat_def = SUMSTATS_DEFINITION[sum_stat].replace("\n", "_").replace("_"," ").title()
        to_df[stat_def] =  [sum_stat, empirical_df[sum_stat], percentile]

    percentile_summary =  pd.DataFrame(to_df,index=["stat_id","empirical","percentile"]).T
    percentile_summary["retained"] = (percentile_summary.percentile.between(2.5,97.5)).astype(int)

    return percentile_summary


def get_max_gap(msa: Msa):
    pattern = re.compile("\w+")
    max_gap = 0
    for seq in msa.get_seq_with_indels():
        new_max = len(max(pattern.split(seq)))
        max_gap = new_max if new_max > max_gap else max_gap
    return max_gap

def fetch_paths(current_path: pathlib.Path):
    if len( n := list(current_path.glob("*.tree"))) == 1:
        tree_path = n[0]

    if len( n := list(current_path.glob("*.fasta"))) == 1:
        msa_path = n[0]

    if tree_path is None or msa_path is None:
        print("no fasta or tree file")
        exit()
    return tree_path, msa_path

def get_retained_data(summary_df: pd.DataFrame):
    return summary_df["retained"].values.tolist()

def get_retained_data_corrected(summary_df: pd.DataFrame):
    return summary_df["retained_with_correction"].values.tolist()



datasets =  ["eggnog","yeast"] # ["adequacy_sims",, "eggnog"]
datasets_str = "".join(datasets)

adequacy_data = {}
retained_summary_statistics_data = {}
retained_summary_statistics_data_corrected = {}

corrected = "corrected"
for current_data in datasets:
    # REMOTE_PATH = REMOTE_PATH.parent.resolve()
    CURRENT_REMOTE_PATH = REMOTE_PATH / current_data
    if not CURRENT_REMOTE_PATH.exists():
        exit(1)

    data_summary = {}
    retained_summary = {}
    retained_summary_corrected = {}

    for dataset_num, prot in enumerate(CURRENT_REMOTE_PATH.iterdir()):
        MAIN_PATH = prot
        TREE_PATH, MSA_PATH = fetch_paths(MAIN_PATH)


        try:
            model_params = raxml_parser.get_substitution_model(MAIN_PATH)
        except (ValueError, IndexError):
            model_params = {
                "inv_prop": -1,
                "gamma_shape": -1
            }

        
        tree = Phylo.read(TREE_PATH, 'newick')
        sum_branch_lengths = tree.total_branch_length()

        remote_emp_msa = Msa(str(MSA_PATH))
        empirical_sum_stats = remote_emp_msa.get_sum_stats()
        empirical_sum_stats = dict(zip(SUMSTATS_LIST, empirical_sum_stats))

        classification = get_classification(MAIN_PATH)
        print(classification)

        max_gap = get_max_gap(remote_emp_msa)

        ADEQUACY_FILES = MAIN_PATH / "adequacy_revised" 
        if not next(ADEQUACY_FILES.glob(f"{classification}_*.csv"), False):
            print(ADEQUACY_FILES)
            continue
        classification_sims = []

        for file_path in ADEQUACY_FILES.glob(f"{classification}_*.csv"):
            classification_sims.append(pd.read_csv(file_path).drop(columns="Unnamed: 0"))
        print(len(classification_sims))
        classification_sims = pd.concat(classification_sims)

        # classification_sims = pd.read_csv(ADEQUACY_FILES / f"{classification}_{corrected}.csv").drop(columns="Unnamed: 0")

        percentile_summary = get_percentile_summary(classification_sims, empirical_sum_stats)
        adequacy_score = sum(percentile_summary.retained)/len(percentile_summary)

        data_summary[prot.stem] = [max_gap, classification, adequacy_score, sum_branch_lengths, model_params["inv_prop"], model_params["gamma_shape"]]
        retained_summary[prot.stem] = [classification, *get_retained_data(percentile_summary)]
        retained_summary_corrected[prot.stem] = [classification, *get_retained_data_corrected(percentile_summary)]

    adequacy_data[current_data] = pd.DataFrame(data_summary)
    adequacy_data[current_data] = adequacy_data[current_data].T.rename(columns={0:"max_gap", 
                                                                                1:"classification",
                                                                                2: "match_percent",
                                                                                3:"sum_branch_lengths",
                                                                                4: "inv_prop",
                                                                                5: "gamma_shape"})
    adequacy_data[current_data] = adequacy_data[current_data].astype({"max_gap": int,
                                                                      "match_percent": float,
                                                                      "sum_branch_lengths":float,
                                                                      "inv_prop": float,
                                                                      "gamma_shape": float})

    retained_summary_statistics_data[current_data] = retained_summary
    retained_summary_statistics_data_corrected[current_data] = retained_summary_corrected

with open(f"assets/adequacy_data_{datasets_str}.pickle", 'wb') as f:
    pickle.dump(file=f, obj=adequacy_data)
with open(f"assets/retained_summary_statistics_data_{datasets_str}.pickle", 'wb') as f:
    pickle.dump(file=f, obj=(retained_summary_statistics_data, retained_summary_statistics_data_corrected))
