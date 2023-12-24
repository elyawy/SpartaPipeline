import pathlib, argparse
from io import StringIO
import pandas as pd
from Bio import Phylo
from tqdm import tqdm
from elyawy.sparta import Simulator, Msa
from indelible_runner import IndelibleCommandline
from raxml_parser import get_substitution_model
from sim_creator import SimConfig



_parser = argparse.ArgumentParser(allow_abbrev=False)
_parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
_parser.add_argument('-o','--output', action='store',metavar="Input folder", type=str, required=True)

_parser.add_argument('-t','--type', action='store',metavar="Type of MSA NT/AA" , type=str, required=True)
_parser.add_argument('-n','--numsim', action='store',metavar="Number of simulations" , type=int, required=True)
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
print(seq_lengths_in_msa)
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

    mas_list = []
    sim_params = sim_config.get_random_sim(num_sims)

    for params in sim_params:
        simulator.init_sim(*params)
        sim_msa = simulator()
        mas_list.append(sim_msa.get_seq_with_indels())
    
    return mas_list
msas_list = init_correction(sim_config, tree_data, NUM_SIMS)

print(msas_list)
def merge_subs_with_indels(msa_list, subsitution_model):
    tree_string = subsitution_model["tree"]
    sparta_organism_order = [i.name for i in Phylo.read(StringIO(tree_string), "newick").get_terminals(order='preorder')]

    fake_msas = []
    for msa in tqdm(msa_list) if VERBOSE else msa_list:
        subsitution_model["length"] = max([len(seq) for seq in msa])
        substitutions = IndelibleCommandline(subsitution_model)
        sub_sim_msa_list = []
        unaligned_msa_list = []
        for idx, sequence in enumerate(msa):
            merged_alignment = ""
            for j,c in enumerate(sequence):
                if c=='-':
                    merged_alignment += "-"
                else:
                    merged_alignment += substitutions[idx][j]
            sub_sim_msa_list.append(merged_alignment)

        sub_sim_msa_list = "".join([f'>{name}\n{seq}\n' for name, seq in zip(sparta_organism_order, sub_sim_msa_list)])
        fake_msas.append(sub_sim_msa_list)

    return fake_msas

fake_msas = merge_subs_with_indels(msas_list, substitution_model)

for idx, msa in enumerate(fake_msas):
    with open(OUTPUT_PATH / f"{idx}.fasta", 'w') as f:
        f.write(msa)