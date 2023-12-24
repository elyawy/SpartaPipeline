# %%
import pathlib,glob
import pandas as pd
# from lightgbm import LGBMRegressor

# %%
REMOTE_PATH = "/groups/pupko/elyawygoda/length_distributions/all_outputs/correction_data"
REMOTE_PATH = pathlib.Path(REMOTE_PATH)
REMOTE_PATH.resolve()

# %%
feature_files  = [f for f in glob.glob(str(REMOTE_PATH / "*_true.parquet.gzip"))]
label_files  = [f for f in glob.glob(str(REMOTE_PATH / "*_realigned.parquet.gzip"))]


# %%
label_col_orig_name = '19'
label_col_new_name = 'label'
feature_dfs = []
for feature_path,label_path in zip(feature_files, label_files):
    print(feature_path)
    feature_df_tmp = pd.read_parquet(feature_path)
    label_df_tmp = pd.read_parquet(label_path)
    feature_df_tmp[label_col_new_name] = label_df_tmp[label_col_orig_name]
    feature_dfs.append(feature_df_tmp)
final_feature_label_df = pd.concat(feature_dfs)


final_feature_label_df.to_parquet("./all_data.parquet.gzip", compression='gzip')

# %

