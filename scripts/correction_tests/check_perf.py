import pathlib, copy, warnings
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn import exceptions
from sklearn.pipeline import Pipeline
warnings.simplefilter("ignore", category=exceptions.ConvergenceWarning)
from lightgbm import plot_importance
from scipy import stats as st
from tqdm import tqdm
from perf_metric import performance_metric, REG_TYPES, get_regression_methods
from elyawy.constants import SUMSTATS_DEFINITION, PARAMS_LIST
# %%
# REMOTE_PATH = "/run/user/1000/gvfs/sftp:host=power.tau.ac.il,user=elyawygoda/groups/pupko/elyawygoda/length_distributions/all_outputs/correction_data"
REMOTE_PATH = pathlib.Path.home() / "Data/correction"
REMOTE_PATH = REMOTE_PATH.resolve()
ALL_DATA_FILE = REMOTE_PATH / "all_correction_data_v2.parquet.gzip"
LABELS = ['label_15', 'label_16', 'label_18', 'label_19']


# %%
if ALL_DATA_FILE.exists():
    final_feature_labels_df = pd.read_parquet(ALL_DATA_FILE)
else:
    print("missing correction data file. exiting.")
    exit(1)
len(final_feature_labels_df)




msa_lengths = []
num_species = []
pearsons_vals = []

print("per dataset cross validation:")
# this is for per dataset cross validation
available_datasets = final_feature_labels_df.dataset.unique().tolist()
label_correlations = {}

np.random.seed(42)

current_datasets = copy.deepcopy(available_datasets)
np.random.shuffle(current_datasets)
# test_size = 100

test_datasets = current_datasets#current_datasets[:test_size]
regressor = get_regression_methods("lgbm")
train_df = final_feature_labels_df
per_label_importance = {}

for LABEL in LABELS:
    per_dataset_metrics = {}
    X_train = train_df.loc[:,:'31']
    stats_columns = [i.replace("\n"," ") for i in SUMSTATS_DEFINITION.values()]
    X_train.columns = list(X_train.columns[:4]) + PARAMS_LIST + stats_columns
    y_train = train_df[LABEL]
    regressor.fit(X_train, y_train)
    per_label_importance[LABEL] = list(regressor.steps[0][1].feature_importances_)


    # plot_importance(regressor.steps[0][1], )
    # # plt.tight_layout(pad=0.5, w_pad=0.3, h_pad=1.0)
    # plt.show()
    continue

    for dataset in tqdm(test_datasets):
        pearson_correlations = []
        test_df = final_feature_labels_df[final_feature_labels_df['dataset'] == dataset]
        test_df = test_df[test_df['indel_model'] == "sim"]
        test_df = test_df[test_df['length_distribution'] == "zipf"]

        # for reg_type in REG_TYPES[1:]:
        # regressor = get_regression_methods(reg_type)


        # if reg_type=="lasso":
        #     X_train = test_df.loc[:,:'31']
        #     y_train = test_df[LABEL]



        X_test = test_df.loc[:,:'31']
        y_test = test_df[LABEL]
        predictions = regressor.predict(X_test)
        pearson_r = st.pearsonr(predictions, y_test)
        pearson_correlations.append(pearson_r[0])

        per_dataset_metrics[dataset] = pearson_correlations

    print(per_dataset_metrics)
    # performance = performance_metric(pearson_correlations)

    # print(performance)
    
    with open(REMOTE_PATH / f"per_dataset_cross_validation_{LABEL}_dif.txt", 'w') as f:
        f.write(str(per_dataset_metrics))

    with open(REMOTE_PATH / f"feature_importace_{LABEL}.txt", 'w') as f:
        f.write(str(per_label_importance))


    print(f"for label {LABEL}")
    print(f"the mean correlation was {np.mean(pearson_correlations)}")
    print(f"the min correlation was {np.min(pearson_correlations)}")
    print(f"the max correlation was {np.max(pearson_correlations)}")

with open(REMOTE_PATH / f"feature_importance.txt", 'w') as f:
    f.write(str(per_label_importance))
