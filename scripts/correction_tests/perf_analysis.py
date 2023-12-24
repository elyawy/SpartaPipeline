import json, pathlib
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from scipy import stats as st
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
# from perf_metric import performance_metric
x = np.array([3, 1, 2])



home_path = pathlib.Path.home()

LABELS = ['label_15', 'label_16', 'label_18', 'label_19']


with open(home_path  / f"Data/correction/feature_importance.txt",'r') as f:
    temp_str =  f.read().replace("'", "\"")
    feature_importance = json.loads(temp_str)

print(feature_importance)

peason_dfs_per_label = {}
bad_datasets_across_labels = []
for label in LABELS:
    print(label)


    with open(home_path  / f"Data/correction/per_dataset_cross_validation_{label}.txt",'r') as f:
        temp_str =  f.read().replace("'", "\"")
        data_dict = json.loads(temp_str)

    data = pd.DataFrame(data_dict).T
    data.columns = ["lasso","lgbm"]

    with open(home_path  / f"Data/correction/per_dataset_cross_validation_{label}_dif.txt",'r') as f:
        temp_str =  f.read().replace("'", "\"")
        data_dict = json.loads(temp_str)

    data["lgbm"] = pd.DataFrame(data_dict).T
    print(data)
    # data["similarity"] = (data.lasso-data.lgbm)**2

    print(st.pearsonr(data["lasso"], data["lgbm"]))
    data["perf_rank"] = pd.cut(data.lasso,10,labels=range(10))

    peason_dfs_per_label[label] = data

    # plt.plot(data[["lasso","lgbm"]].sort_values(by="lgbm"), label=["lasso","lgbm"])
    # plt.legend()
    # plt.show()

#     bad_data = data[data.lasso > 0.85]
#     bad_datasets_across_labels.append(set(bad_data.index))

# bad_datasets_intersection = set.intersection(*bad_datasets_across_labels)
# print(len(bad_datasets_intersection))

REMOTE_PATH = pathlib.Path.home() / "Data/correction"
REMOTE_PATH = REMOTE_PATH.resolve()
ALL_DATA_FILE = REMOTE_PATH / "all_correction_data_v2.parquet.gzip"

if ALL_DATA_FILE.exists():
    final_feature_labels_df = pd.read_parquet(ALL_DATA_FILE)
else:
    print("missing correction data file. exiting.")
    exit(1)
len(final_feature_labels_df)

print(final_feature_labels_df.columns)
for label in LABELS:
    rank_statistics = []
    data = peason_dfs_per_label[label]

    rank_groups = data.groupby("perf_rank").groups
    for key,values in data.groupby("perf_rank").groups.items():
        group_label_stats = []
        print(key)
        for dataset in values:
            group_feature_labels_df = final_feature_labels_df[final_feature_labels_df.dataset == dataset]
            group_feature_labels_df = group_feature_labels_df[group_feature_labels_df.length_distribution == "zipf"]
            group_feature_labels_df = group_feature_labels_df[group_feature_labels_df.indel_model == "sim"]

            group_feature_labels_df = group_feature_labels_df.select_dtypes('number')
            stats_series = pd.DataFrame(data=[(key, dataset)],
                                        columns=["perf_rank","dataset"])
            feature_labels_df_std = group_feature_labels_df.std()
            stats_series[feature_labels_df_std.index] = feature_labels_df_std
            group_label_stats.append(stats_series)
            # group_label_stats.append(group_feature_labels_df.value_counts())
        if not values.empty:
            rank_statistics.append(pd.concat(group_label_stats, axis=0))

    rank_statistics = pd.concat(rank_statistics).reset_index(drop=True)

    var_corr = rank_statistics.loc[:,"dataset":].corrwith(rank_statistics["perf_rank"])

    print(var_corr[var_corr > 0.5])
    print(var_corr[var_corr < -0.5])



# st.pearsonr(rank_statistics["perf_rank"], rank_statistics["5"])
# plt.scatter(rank_statistics["perf_rank"], rank_statistics["5"])
# plt.show()

# st.pearsonr(rank_statistics["perf_rank"], rank_statistics["label_13"])
# plt.scatter(rank_statistics["perf_rank"], rank_statistics["label_13"])
# plt.show()
# clf = LinearRegression()
# scaler = StandardScaler()
# scaled_rank_stats = scaler.fit_transform(rank_statistics.loc[:,:"count_mean"])
# clf.fit(scaled_rank_stats, rank_statistics.perf_rank)

# print(st.pearsonr(rank_statistics.perf_rank, clf.predict(scaled_rank_stats)))
# plt.scatter(rank_statistics.perf_rank, clf.predict(scaled_rank_stats))
# plt.show()

# bad_datasets_df = final_feature_labels_df[final_feature_labels_df["dataset"].isin(bad_datasets_intersection)]
# print(final_feature_labels_df[final_feature_labels_df["dataset"].isin(bad_datasets_intersection)])

# bad_datasets_df = bad_datasets_df[bad_datasets_df.length_distribution == "zipf"]
# bad_datasets_df = bad_datasets_df[bad_datasets_df.indel_model == "sim"]
# print(bad_datasets_df)
# for dataset in bad_datasets_intersection:
#     print(dataset)
#     print(bad_datasets_df[bad_datasets_df.dataset == dataset]["label_16"].value_counts())

