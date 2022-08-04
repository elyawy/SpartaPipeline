import os, pathlib, tempfile, pickle, argparse, warnings, copy
from io import StringIO
import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Align.Applications import MafftCommandline
from indelible_runner import IndelibleCommandline
from raxml_parser import get_substitution_model
from sklearn import linear_model, model_selection, exceptions
warnings.simplefilter("ignore", category=exceptions.ConvergenceWarning)
from scipy.stats import pearsonr
from sim_creator import SimConfig
import Sparta

gr_parser = argparse.ArgumentParser(allow_abbrev=False)
gr_parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
# gr_parser.add_argument('-c','--config', action='store',metavar="Simulation config" , type=str, required=True)
gr_parser.add_argument('-t','--type', action='store',metavar="Type of MSA NT/AA" , type=str, required=True)
gr_parser.add_argument('-n','--numsim', action='store',metavar="Number of simulations" , type=int, required=True)
# gr_parser.add_argument('-s','--seed', action='store',metavar="Simulation config" , type=int, required=False)
gr_parser.add_argument('-l','--lengthdist', action='store',metavar="Simulation config" , type=str, required=True)
gr_parser.add_argument('-m','--model', action='store',metavar="Simulation config" , type=str, required=True)


args = gr_parser.parse_args()

MAIN_PATH = args.input
MODE = args.type
NUM_SIMS = args.numsim
LENGTH_DISTRIBUTION = args.lengthdist
INDEL_MODEL = args.model

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

full_correction_path = pathlib.Path(MAIN_PATH, "correction")
try:
    os.mkdir(full_correction_path)
except:
    print("correction folder exists already")




true_msa = Sparta.Msa(MSA_PATH)
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
    simulator = Sparta.Sim(tree_data["path"])

    correction_list = []
    correction_list_sum_stats = []

    sim_params_correction = sim_config.get_random_sim(num_sims)

    for params in sim_params_correction:
        numeric_params = [params[0],params[2][0], params[3][0], params[4], params[5]]
        simulator.init_sim(*params)

        sim_msa = simulator.run_sim()
        correction_list.append(sim_msa.get_seq_with_indels())
        correction_list_sum_stats.append(numeric_params + sim_msa.get_sum_stats())
    
    return correction_list, correction_list_sum_stats
msas_list, sum_stats_list = init_correction(sim_config, tree_data, args.numsim)

def get_regressors(msa_list, sum_stats_list, subsitution_model):
    tree_string = subsitution_model["tree"]
    sparta_organism_order = [i.name for i in Phylo.read(StringIO(tree_string), "newick").get_terminals(order='preorder')]

    realigned_sum_stats = []
    for msa in msa_list:
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
            mafft_cline = MafftCommandline(input=tempf.name)
            realigned_msa, stderr = mafft_cline()

        realigned_msa = [s[s.index("\n"):].replace("\n","") for s in realigned_msa.split(">")[1:]]
        realigned_sum_stats.append(Sparta.Msa(realigned_msa).get_sum_stats())

    def compute_regressors(true_stats, corrected_stats):
        X = np.array(true_stats, dtype=float)
        Y = np.array(corrected_stats, dtype=float).T
        X_train = (X-X.mean())/X.std()
        Y_train = (Y-Y.mean())/Y.std()
        reg = linear_model.Lasso()
        parameters = {'alpha':np.logspace(-7,4,20)}
        clf_lassocv = model_selection.GridSearchCV(estimator = reg,
                                    param_grid = parameters, cv=3,
                                    scoring = 'neg_mean_squared_error')
        regressors = []
        performance_metrics = []
        for y in Y_train:
            clf_lassocv.fit(X_train, y)
            saved_estimator = copy.deepcopy(clf_lassocv)
            regressors.append(saved_estimator)
            
            Y_pred = clf_lassocv.predict(X_train)
            r_val, p_val = pearsonr(Y_pred,y)
            performance_metrics.append({
                'pearsonr': r_val,
                'p_val': p_val,
                'mean_test_score': np.min(np.sqrt(-clf_lassocv.cv_results_['mean_test_score']))
            })
        return regressors, performance_metrics

    regressors, performance = compute_regressors(sum_stats_list, realigned_sum_stats)        

    return regressors, performance
regressors, performance = get_regressors(msas_list, sum_stats_list, substitution_model)

if MAIN_PATH is None:
    for weight in [reg.best_model_.coef_ for reg in regressors]:
        print(",".join([str(w) for w in weight]))
else:
    full_correction_path = os.path.join(MAIN_PATH, "correction")
    path_joiner = lambda x: os.path.join(full_correction_path, x)
    try:
        os.mkdir(full_correction_path)
    except:
        print("correction folder exists already")
    pickle.dump(regressors, open(path_joiner(f'regressors_{LENGTH_DISTRIBUTION}_{INDEL_MODEL}'), 'wb'))
    pd.DataFrame(performance).to_csv(path_joiner(f'regression_performance_{LENGTH_DISTRIBUTION}_{INDEL_MODEL}.csv'))


    




