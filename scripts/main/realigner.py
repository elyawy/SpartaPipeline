import pathlib, tempfile, warnings, copy
import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin# define the transformer
from sklearn import linear_model, model_selection, exceptions
from sklearn.pipeline import Pipeline
warnings.simplefilter("ignore", category=exceptions.ConvergenceWarning)
from scipy.stats import pearsonr

from io import StringIO
from Bio import Phylo
from Bio.Align.Applications import MafftCommandline
from elyawy.sparta import Msa
from indelible_runner import IndelibleCommandline
from raxml_parser import get_substitution_model

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


def realign_sim_msa(msa: Msa, tree_string: str, alphabet: str, model_path: pathlib.Path=""):




    sparta_organism_order = [i.name for i in Phylo.read(StringIO(tree_string), "newick").get_terminals(order='preorder')]


    indels = np.array([list(x) for x in msa.get_seq_with_indels()])
    msa_len = len(indels[0])

    # prepare indelible control file for subsitutions:
    substitution_model = get_substitution_model(str(model_path)) if alphabet == "NT" else {}
    substitution_model["tree"] = tree_string
    substitution_model["mode"] = "nuc" if alphabet == "NT" else "amino"
    substitution_model["length"] = msa_len
    substitutions = IndelibleCommandline(substitution_model)
    
    new_msa = np.ascontiguousarray([list(x.strip()) for x in substitutions])
    np.putmask(new_msa, indels == "-", "-")
    
    new_msa = new_msa.view(f'U{msa_len}').ravel()
    new_msa = "".join([f'>{name}\n{seq}\n' for name, seq in zip(sparta_organism_order, new_msa)])
    new_msa = new_msa.encode()

    with tempfile.NamedTemporaryFile(suffix='.fasta') as tempf:
        tempf.write(new_msa)
        tempf.seek(0)
        mafft_cline = MafftCommandline(input=tempf.name)
        realigned_msa, stderr = mafft_cline()

    realigned_msa = [s[s.index("\n"):].replace("\n","") for s in realigned_msa.split(">")[1:]]

    return Msa(realigned_msa)


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
