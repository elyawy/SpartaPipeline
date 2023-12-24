import pathlib, pickle, uuid
import numpy as np
from elyawy.sparta import Simulator
import realigner as rl
# EXPERIMENT 3 no ML

WORKING_DIR = pathlib.Path.home() / "Data/adequacy_experiments"
pathlib.Path.home()
WORKING_DIR.exists()
CURRENT_EXPERIMENT_PATH = WORKING_DIR /  "experiment_3"

TREE_PATH = None
if len( n := list(CURRENT_EXPERIMENT_PATH.glob("*.tree"))) == 1:
    TREE_PATH = n[0]
TREE_PATH
sim = Simulator(TREE_PATH)

NUMBER_OF_SAMPLES = 1000
sim_msas_stats = np.zeros(shape=(NUMBER_OF_SAMPLES,27))
tree_str = TREE_PATH.read_text()

empirical_params = np.array([[470, 0.03811097603185668, 0.03811097603185668, 'zipf',
                              1.3349875044949302, 1.3349875044949302]])



for i in range(NUMBER_OF_SAMPLES):
    sim_msa = sim()
    sim_msa_realigned = rl.realign_sim_msa(sim_msa, tree_str, "AA")
    sim_msas_stats[i,:] = sim_msa_realigned.get_sum_stats()

    print(i)

with open(CURRENT_EXPERIMENT_PATH / f"{uuid.uuid1()}.pickle", 'wb') as f:
    pickle.dump(sim_msas_stats, f)
