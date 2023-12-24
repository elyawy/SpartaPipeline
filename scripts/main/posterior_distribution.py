import pathlib, argparse
import pandas as pd
from elyawy.sparta import Simulator
import realigner as rl




_parser = argparse.ArgumentParser(allow_abbrev=False)
_parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
# _parser.add_argument('-c','--config', action='store',metavar="Simulation config" , type=str, required=True)
_parser.add_argument('-n','--numsim', action='store',metavar="Number of simulations" , type=int, required=True)
# _parser.add_argument('-s','--seed', action='store',metavar="Simulation config" , type=int, required=False)
_parser.add_argument('-t','--type', action='store',metavar="Type of MSA NT/AA" , type=str, required=True)

_parser.add_argument('-nc','--no-correction', action='store_false') # default is with correction


args = _parser.parse_args()

WORKING_PATH = pathlib.Path(args.input)
NUM_SIMS = 10000
TREE_PATH = next(WORKING_PATH.glob("*.tree"))
CORRECTION = args.no_correction # default is with correction
MSA_TYPE = args.type






sim = Simulator(TREE_PATH)




parameters_file = next(WORKING_PATH.glob("*.parquet"))


posterior_parameters = pd.read_parquet(parameters_file)

sampled_parameters = posterior_parameters.sample(NUM_SIMS, replace=True).reset_index(drop=True)


# steps:
#generate simulations
tree_str = TREE_PATH.read_text()
for i,params in enumerate(sampled_parameters):
        sim.init_sim(*params)
        msa = sim()
        if CORRECTION:
            msa = rl.realign_sim_msa(msa, tree_str, MSA_TYPE)
        sim_msas_stats[i,:] = msa.get_sum_stats()
