# %%
import pathlib
import pandas as pd
from lightgbm import LGBMRegressor
import numpy as np
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split, GridSearchCV
from scipy import stats as st


# %%
# REMOTE_PATH = "/run/user/1000/gvfs/sftp:host=power.tau.ac.il,user=elyawygoda/groups/pupko/elyawygoda/length_distributions/all_outputs/correction_data"
REMOTE_PATH = "/groups/pupko/elyawygoda/length_distributions/all_outputs/correction_data"
REMOTE_PATH = pathlib.Path(REMOTE_PATH)
REMOTE_PATH.resolve()
ALL_DATA_FILE = pathlib.Path("elyawygoda/groups/pupko/elyawygoda/length_distributions/all_outputs/all_correction_data.parquet.gzip")
LABEL_COLUMNS = [f"label_{i}" for i in range(0,27)]

# %%
def fetch_all_data_df_from_server() -> pd.DataFrame:
    feature_files  = sorted(list(REMOTE_PATH.glob("*_true.parquet.gzip")))
    label_files  = sorted(list(REMOTE_PATH.glob("*_realigned.parquet.gzip")))
    for feature_path,label_path in zip(feature_files, label_files):
        dataset_a = "_".join(feature_path.stem.split("_")[:4])
        dataset_b = "_".join(label_path.stem.split("_")[:4])
        assert dataset_a == dataset_b
        
    

    feature_dfs = []
    for feature_path,label_path in zip(feature_files, label_files):
        feature_df_tmp = pd.read_parquet(feature_path)
        label_df_tmp = pd.read_parquet(label_path)
        feature_df_tmp[LABEL_COLUMNS] = label_df_tmp
        dataset, lendist, indel_model, _ = feature_path.stem.rsplit("_", maxsplit=3)
        feature_df_tmp[["dataset", "length_distribution", "indel_model"]] = dataset, lendist, indel_model
        feature_dfs.append(feature_df_tmp)
    final_feature_labels_df = pd.concat(feature_dfs)

    return final_feature_labels_df


# %%
if ALL_DATA_FILE.exists():
    final_feature_labels_df = pd.read_parquet(ALL_DATA_FILE)
else:
    final_feature_labels_df = fetch_all_data_df_from_server()
    final_feature_labels_df.to_parquet(ALL_DATA_FILE, compression="gzip")
len(final_feature_labels_df)
# %%
bad_labels = ['label_15', 'label_16', 'label_18', 'label_19']
# %%
regressors = {}
for label in bad_labels:
    labels = final_feature_labels_df[label]
    features = final_feature_labels_df.loc[:,:'31']
    X_train, X_test, y_train, y_test = train_test_split(features,
                                                        labels,
                                                        train_size=0.9)
    regressor = GridSearchCV(LGBMRegressor(), 
                             param_grid={
                              'boosting_type': ['gbdt', 'dart', 'rf'],
                              'num_leaves': [32,40],
                              'n_estimators': [100, 200, 300],
                              'reg_alpha': np.logspace(-3,4,5),
                              'reg_lambda': np.logspace(-3,4,5),
                              'random_state': [42]
                             },
                             scoring='neg_mean_squared_error',
                            )
    regressor.fit(X_train, y_train)
    regressors[label] = regressor.best_estimator_
    predictions = regressor.predict(X_test)
    pearson_r = st.pearsonr(predictions,y_test )
    print(f"correlation for {label} is: {pearson_r[0]:.3f}")




for reg in regressors.values():
    with open("/groups/pupko/elyawygoda/length_distributions/all_outputs/best_estimators_bad_labels.txt", 'w') as f:
        f.write(str(reg.get_params()))
