import pathlib, tempfile, argparse, sys
from io import StringIO
import pandas as pd
from Bio import Phylo
from Bio.Align.Applications import MafftCommandline
from tqdm import tqdm
from elyawy.sparta import Simulator, Msa
from indelible_runner import IndelibleCommandline
from raxml_parser import get_substitution_model
from sim_creator import SimConfig



_parser = argparse.ArgumentParser(allow_abbrev=False)
_parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
_parser.add_argument('-o','--output', action='store',metavar="Input folder", type=str, required=True)

# _parser.add_argument('-c','--config', action='store',metavar="Simulation config" , type=str, required=True)
_parser.add_argument('-t','--type', action='store',metavar="Type of MSA NT/AA" , type=str, required=True)
_parser.add_argument('-n','--numsim', action='store',metavar="Number of simulations" , type=int, required=True)
# _parser.add_argument('-s','--seed', action='store',metavar="Simulation config" , type=int, required=False)
_parser.add_argument('-l','--lengthdist', action='store',metavar="Simulation config" , type=str, required=True)
_parser.add_argument('-m','--model', action='store',metavar="Simulation config" , type=str, required=True)
_parser.add_argument('-v','--verbose', action='store_true')


args = _parser.parse_args()

MAIN_PATH = pathlib.Path(args.input).resolve()
OUTPUT_PATH = pathlib.Path(args.output).resolve()
MODE = args.type
NUM_SIMS = args.numsim
LENGTH_DISTRIBUTION = args.lengthdist
INDEL_MODEL = args.model
VERBOSE = args.verbose

def fetch_paths(current_path: pathlib.Path):
    if len( n := list(current_path.glob("*.tree"))) == 1:
        tree_path = n[0]

    if len( n := list(current_path.glob("*.fasta"))) == 1:
        msa_path = n[0]

    if tree_path is None or msa_path is None:
        print("no fasta or tree file")
        exit()
    return tree_path, msa_path
TREE_PATH, MSA_PATH = fetch_paths(MAIN_PATH)
TREE_PATH, MSA_PATH = str(TREE_PATH), str(MSA_PATH)


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
substitution_model = get_substitution_model(str(MAIN_PATH)) if MODE == "NT" else {}
substitution_model["tree"] = tree_string
substitution_model["mode"] = "nuc" if MODE == "NT" else "amino"

def init_correction(sim_config, tree_data, num_sims):
    simulator = Simulator(tree_data["path"])

    correction_list = []
    correction_list_sum_stats = []

    sim_params_correction = sim_config.get_random_sim(num_sims)

    for params in sim_params_correction:
        numeric_params = [params[0],params[1], params[2], params[4], params[5]]
        simulator.init_sim(*params)
        sim_msa = simulator()

        correction_list.append(sim_msa.get_seq_with_indels())
        correction_list_sum_stats.append(numeric_params + sim_msa.get_sum_stats())
    
    return correction_list, correction_list_sum_stats
msas_list, sum_stats_list = init_correction(sim_config, tree_data, NUM_SIMS)


def get_realigned_stats(msa_list, subsitution_model):
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
            mafft_cline = MafftCommandline(input=tempf.name)
            realigned_msa, stderr = mafft_cline()

        realigned_msa = [s[s.index("\n"):].replace("\n","") for s in realigned_msa.split(">")[1:]]
        realigned_sum_stats.append(Msa(realigned_msa).get_sum_stats())


    return realigned_sum_stats
realigned_sum_stats = get_realigned_stats(msas_list, substitution_model)


saving_str = f"{MAIN_PATH.stem}_{LENGTH_DISTRIBUTION}_{INDEL_MODEL}"

print("saving stats...")
true_stats = pd.DataFrame(sum_stats_list)
true_stats.columns = map(str, range(len(true_stats.columns)))
true_stats.to_parquet(OUTPUT_PATH / f"{saving_str}_true.parquet.gzip", compression='gzip', index=False)

realigned_stats = pd.DataFrame(realigned_sum_stats)
realigned_stats.columns = map(str, realigned_stats.columns)
realigned_stats.to_parquet(OUTPUT_PATH / f"{saving_str}_realigned.parquet.gzip", compression='gzip', index=False)