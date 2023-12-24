import os, pathlib, tempfile, pickle, argparse, warnings, copy
from io import StringIO
import numpy as np
import pandas as pd
from Bio import Phylo
from tqdm import tqdm
from indelible_runner import IndelibleCommandline
from raxml_parser import get_substitution_model
from sklearn.base import BaseEstimator, TransformerMixin# define the transformer
from sklearn import linear_model, model_selection, exceptions
from sklearn.pipeline import Pipeline
warnings.simplefilter("ignore", category=exceptions.ConvergenceWarning)
from scipy.stats import pearsonr
from sim_creator import SimConfig
from elyawy.sparta import Simulator, Msa
from aligner_interface import Aligner

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
# _parser.add_argument('-c','--config', action='store',metavar="Simulation config" , type=str, required=True)
_parser.add_argument('-t','--type', action='store',metavar="Type of MSA NT/AA" , type=str, required=True)
_parser.add_argument('-n','--numsim', action='store',metavar="Number of simulations" , type=int, required=True)
# _parser.add_argument('-s','--seed', action='store',metavar="Simulation config" , type=int, required=False)
_parser.add_argument('-a','--aligner', action='store',metavar="Alignment program to use" , type=str, required=True)

_parser.add_argument('-l','--lengthdist', action='store',metavar="Simulation config" , type=str, required=True)
_parser.add_argument('-m','--model', action='store',metavar="Simulation config" , type=str, required=True)
_parser.add_argument('-k','--keep-stats', action='store_true')
_parser.add_argument('-v','--verbose', action='store_true')


args = _parser.parse_args()

MAIN_PATH = args.input
MODE = args.type
NUM_SIMS = args.numsim
ALIGNER = Aligner(args.aligner)
LENGTH_DISTRIBUTION = args.lengthdist
INDEL_MODEL = args.model
KEEP_STATS = args.keep_stats
VERBOSE = args.verbose


TREE_PATH = None
MSA_PATH = None
for file in os.listdir(MAIN_PATH):
    if ".tree" in file:
        TREE_PATH = os.path.join(MAIN_PATH, file)
    if ".fasta" in file:
        MSA_PATH = os.path.join(MAIN_PATH, file)

if TREE_PATH is None or MSA_PATH is None:
    print("no fasta or tree file")
    exit()


true_msa = Msa(MSA_PATH)
seq_lengths_in_msa = [true_msa.get_shortest_seq_length(), true_msa.get_longest_seq_length()]

sim_config = SimConfig(seq_lengths=seq_lengths_in_msa,
                       len_dist=LENGTH_DISTRIBUTION,
                       indel_model=INDEL_MODEL)
# np.random.seed(int(sim_config.seed))

with open(TREE_PATH) as treef:
        tree_string = treef.read()

tree_data = {
    "string": tree_string,
    "path": TREE_PATH
}

# prepare indelible control file for subsitutions:
substitution_model = get_substitution_model(MAIN_PATH) if MODE == "NT" else {}
substitution_model["tree"] = tree_string
substitution_model["mode"] = "nuc" if MODE == "NT" else "amino"

def init_correction(sim_config, tree_data, num_sims):
    simulator = Simulator(tree_data["path"])

    correction_list = []
    correction_list_sum_stats = []

    sim_params_correction = sim_config.get_random_sim(num_sims)

    for idx,params in tqdm(enumerate(sim_params_correction)) if VERBOSE else sim_params_correction:
        sim_stats = [0]*27
        while not any(sim_stats):
            numeric_params = [params[0],params[1], params[2], params[4], params[5]]
            simulator.init_sim(*params)
            sim_msa = simulator()
            sim_stats = sim_msa.get_sum_stats()
            if not any(sim_stats):
                replacement_params = sim_config.get_random_sim(1)
                sim_params_correction[idx] = replacement_params[0]
        correction_list.append(sim_msa.get_seq_with_indels())
        correction_list_sum_stats.append(numeric_params + sim_stats)    

    return correction_list, correction_list_sum_stats
msas_list, sum_stats_list = init_correction(sim_config, tree_data, NUM_SIMS)


def get_regressors(msa_list, sum_stats_list, subsitution_model):
    tree_string = subsitution_model["tree"]
    sparta_organism_order = [i.name for i in Phylo.read(StringIO(tree_string), "newick").get_terminals(order='preorder')]

    realigned_sum_stats = []
    for msa in tqdm(msa_list) if VERBOSE else msa_list:
        subsitution_model["length"] = max([len(seq) for seq in msa])
        substitutions = IndelibleCommandline(subsitution_model)
        # sub_sim_msa_list = []
        unaligned_msa_list = []
        for idx, sequence in enumerate(msa):
            merged_alignment = ""
            for j,c in enumerate(sequence):
                if c=='-':
                    merged_alignment += "-"
                else:
                    merged_alignment += substitutions[idx][j]
            # sub_sim_msa_list.append(merged_alignment)
            unaligned_msa_list.append(merged_alignment.replace('-',''))

        sim_fasta_unaligned = "".join([f'>{name}\n{seq}\n' for name, seq in zip(sparta_organism_order, unaligned_msa_list)])

        sim_fasta_unaligned = sim_fasta_unaligned.encode()

        with tempfile.NamedTemporaryFile(suffix='.fasta') as tempf:
            tempf.write(sim_fasta_unaligned)
            tempf.seek(0)
            ALIGNER.set_input_file(tempf.name, tree_file=TREE_PATH)
            realigned_msa = ALIGNER.get_realigned_msa()
        
        print(realigned_msa)
        realigned_msa = [s[s.index("\n"):].replace("\n","") for s in realigned_msa.split(">")[1:]]
        realigned_sum_stats.append(Msa(realigned_msa).get_sum_stats())

    def compute_regressors(true_stats, corrected_stats):
        X = np.array(true_stats, dtype=float)
        Y = np.array(corrected_stats, dtype=float).T

        reg = linear_model.Lasso()
        parameters = {'alpha':np.logspace(-7,4,20)}
        clf_lassocv = model_selection.GridSearchCV(estimator = reg,
                                    param_grid = parameters, cv=3,
                                    scoring = 'neg_mean_squared_error')
        regression_pipline = Pipeline([("scaler", StandardMemoryScaler()),('regression', clf_lassocv)])
        regressors = []
        performance_metrics = []
        for y in Y:
            regression_pipline.fit(X, y)
            saved_estimator = copy.deepcopy(regression_pipline)
            regressors.append(saved_estimator)
            
            Y_pred = regression_pipline.predict(X)
            r_val, p_val = pearsonr(Y_pred,y)
            performance_metrics.append({
                'pearsonr': r_val,
                'p_val': p_val,
                'mean_test_score': np.min(np.sqrt(-clf_lassocv.cv_results_['mean_test_score']))
            })
        return regressors, performance_metrics

    regressors, performance = compute_regressors(sum_stats_list, realigned_sum_stats)        

    return regressors, performance, realigned_sum_stats
regressors, performance, realigned_sum_stats = get_regressors(msas_list, sum_stats_list, substitution_model)


if MAIN_PATH is None:
    for weight in [reg.best_model_.coef_ for reg in regressors]:
        print(",".join([str(w) for w in weight]))
else:
    full_correction_path = os.path.join(MAIN_PATH, f"{args.aligner}_correction")
    path_joiner = lambda x: os.path.join(full_correction_path, x)
    try:
        os.mkdir(full_correction_path)
    except:
        print("correction folder exists already")
    pickle.dump(regressors, open(path_joiner(f'regressors_{LENGTH_DISTRIBUTION}_{INDEL_MODEL}'), 'wb'))
    pd.DataFrame(performance).to_csv(path_joiner(f'regression_performance_{LENGTH_DISTRIBUTION}_{INDEL_MODEL}.csv'))

if KEEP_STATS:
    print("saving stats...")
    true_stats = pd.DataFrame(sum_stats_list)
    true_stats.columns = map(str, range(len(true_stats.columns)))
    true_stats.to_parquet("true_stats.parquet.gzip", compression='gzip', index=False)

    realigned_stats = pd.DataFrame(realigned_sum_stats)
    realigned_stats.columns = map(str, realigned_stats.columns)
    realigned_stats.to_parquet("realigned_stats.parquet.gzip", compression='gzip', index=False)

    infered_stats = np.array([regressor.predict(true_stats.values).T for regressor in regressors])
    infered_stats = pd.DataFrame(infered_stats.T, columns=map(str, range(27)))
    infered_stats.to_parquet("infered_stats.parquet.gzip", compression='gzip', index=False)



