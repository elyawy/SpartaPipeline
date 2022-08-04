import pickle
path_to_file = "/home/elyawy/Development/Msc/sandbox/seqal_chrI_87389-87501_+_YAL030W.aln/correction/regressors_geometric"
with open(path_to_file, 'rb') as f:
    val = pickle.load(f)

for v in val:
    print(id(v))
