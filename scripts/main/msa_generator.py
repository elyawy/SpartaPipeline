import pathlib, argparse, tempfile, shutil
from io import StringIO

import numpy as np

from elyawy.sparta import Simulator, Msa
from elyawy.constants import length_distributions, indel_models
from sim_creator import SimConfig
from raxml_parser import get_substitution_model
from indelible_runner import IndelibleCommandline

from Bio.Align.Applications import MafftCommandline
from Bio import Phylo



_parser = argparse.ArgumentParser(allow_abbrev=False)

_parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
_parser.add_argument('-o','--target', action='store',metavar="Target folder", type=str, required=True)
_parser.add_argument('-t','--type', action='store',metavar="Type of MSA NT/AA" , type=str, required=True)
_parser.add_argument('-n','--numsim', action='store',metavar="Number of simulations" , type=int, required=True)
_parser.add_argument('-l','--lengthdist', action='store',metavar="Length distribution" , type=str, required=False)
_parser.add_argument('-m','--model', action='store',metavar="indel model" , type=str, required=False)



args = _parser.parse_args()

NUM_SIMS = args.numsim
LENGTH_DISTRIBUTION = args.lengthdist
if not args.lengthdist:
    LENGTH_DISTRIBUTION = np.random.choice(length_distributions, 1)[0]
INDEL_MODEL = args.model
if not args.model:
    INDEL_MODEL = np.random.choice(indel_models, 1)[0]

MODE = args.type

INPUT_PATH = pathlib.Path(args.input).resolve()

TARGET_PATH = pathlib.Path(args.target).resolve() / f"{INPUT_PATH.stem}_{LENGTH_DISTRIBUTION}_{INDEL_MODEL}"

TARGET_PATH.mkdir(exist_ok=True)
TREE_PATH = None
MSA_PATH = None
for file in INPUT_PATH.iterdir():
    if ".tree" in str(file):
        TREE_PATH = str(INPUT_PATH / file)
    if ".fasta" in str(file):
        MSA_PATH = str(INPUT_PATH / file)

if TREE_PATH is None or MSA_PATH is None:
    print("no fasta or tree file")
    exit()


true_msa = Msa(MSA_PATH)
sim = Simulator(TREE_PATH)


true_msa_sum_stats = np.array(true_msa.get_sum_stats())
seq_lengths_in_msa = [true_msa.get_shortest_seq_length(), true_msa.get_longest_seq_length()]


sim_config = SimConfig(len_dist=LENGTH_DISTRIBUTION,
                       seq_lengths=seq_lengths_in_msa,
                       indel_model=INDEL_MODEL)

# SEED = int(sim_config.seed)
# np.random.seed(SEED)
# sim.set_seed(SEED)

sim_params = sim_config.get_random_sim(NUM_SIMS)

with open(TREE_PATH) as treef:
        tree_string = treef.read()

tree_data = {
    "string": tree_string,
    "path": TREE_PATH
}

substitution_model = get_substitution_model(INPUT_PATH) if MODE == "NT" else {}
substitution_model["tree"] = tree_string
substitution_model["mode"] = "nuc" if MODE == "NT" else "amino"

msa_list = []
for params in sim_params:

    numeric_params = [params[0],params[1], params[2], params[4], params[5]]
    sim.init_sim(*params)
    sim_msa = sim()
    msa_list.append(sim_msa.get_seq_with_indels())

tree_string = substitution_model["tree"]
sparta_organism_order = [i.name for i in Phylo.read(StringIO(tree_string), "newick").get_terminals(order='preorder')]


for msa_num, msa in enumerate(msa_list):
    substitution_model["length"] = max([len(seq) for seq in msa])
    substitutions = IndelibleCommandline(substitution_model)
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
        unaligned_msa_list.append(merged_alignment.replace('-',''))

    sim_fasta_aligned = "".join([f'>{name}\n{seq}\n' for name, seq in zip(sparta_organism_order, sub_sim_msa_list)])

    with open(TARGET_PATH / f"true_msa_{msa_num}.txt", 'w') as f:
        f.write(sim_fasta_aligned)

    sim_fasta_unaligned = "".join([f'>{name}\n{seq}\n' for name, seq in zip(sparta_organism_order, unaligned_msa_list)])
    sim_fasta_unaligned = sim_fasta_unaligned.encode()
    with tempfile.NamedTemporaryFile(suffix='.fasta') as tempf:
        tempf.write(sim_fasta_unaligned)
        tempf.seek(0)
        mafft_cline = MafftCommandline(input=tempf.name)
        realigned_msa, stderr = mafft_cline()
    
    with open(TARGET_PATH / f"mafft_realigned_{msa_num}.fasta", 'w') as f:
        f.write(realigned_msa)

    shutil.copy2(TREE_PATH, TARGET_PATH)