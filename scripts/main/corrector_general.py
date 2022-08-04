import os, tempfile, argparse, pathlib
from io import StringIO
import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Align.Applications import MafftCommandline
from Bio.Align.Applications import PrankCommandline

from indelible_runner import IndelibleCommandline
from raxml_parser import get_substitution_model
from sim_creator import SimConfig
from elyawy.constants import PARAMS_LIST, SUMSTATS_DEFINITION
import Sparta

gr_parser = argparse.ArgumentParser(allow_abbrev=False)
gr_parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
# gr_parser.add_argument('-c','--config', action='store',metavar="Simulation config" , type=str, required=True)
gr_parser.add_argument('-t','--type', action='store',metavar="Type of MSA NT/AA" , type=str, required=True)
gr_parser.add_argument('-n','--numsim', action='store',metavar="Number of simulations" , type=int, required=True)
# gr_parser.add_argument('-s','--seed', action='store',metavar="Simulation config" , type=int, required=False)
gr_parser.add_argument('-l','--lengthdist', action='store',metavar="Simulation config" , type=str, required=True)
gr_parser.add_argument('-m','--model', action='store',metavar="Simulation config" , type=str, required=True)
gr_parser.add_argument('-a','--alignmentprogram', action='store',metavar="Alignment program" , type=str, required=True)


args = gr_parser.parse_args()

MAIN_PATH = args.input
MODE = args.type
NUM_SIMS = args.numsim
LENGTH_DISTRIBUTION = args.lengthdist
INDEL_MODEL = args.model
ALIGNMENT_PROGRAM = args.alignmentprogram

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

# generate #num_sims msas without substitutions according to sim_config and tree_data
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
msas_list, sum_stats_list = init_correction(sim_config, tree_data, NUM_SIMS)

tree_string = substitution_model["tree"]
sparta_organism_order = [i.name for i in Phylo.read(StringIO(tree_string), "newick").get_terminals(order='preorder')]

def add_substitutions(msa, true_msas_path):
    substitution_model["length"] = max([len(seq) for seq in msa])
    substitutions = IndelibleCommandline(substitution_model)
    
    sub_sim_msa_list = []
    for idx, sequence in enumerate(msa):
        merged_alignment = ""
        for j,c in enumerate(sequence):
            if c=='-':
                merged_alignment += "-"
            else:
                merged_alignment += substitutions[idx][j]
        sub_sim_msa_list.append(merged_alignment)

    sim_fasta_aligned = "".join([f'>{name}\n{seq}\n' for name, seq in zip(sparta_organism_order, sub_sim_msa_list)])
    
    run_id = f"{INDEL_MODEL}_{LENGTH_DISTRIBUTION}_{ALIGNMENT_PROGRAM}"
    file_paths = filter(lambda x: run_id in x.stem ,true_msas_path.iterdir())
    file_index_list = [int(file_path.stem.split("_")[-1]) for file_path in file_paths]
    new_file_index = 0 if not file_index_list else max(file_index_list) + 1
    print(new_file_index)
    file_id = f"{INDEL_MODEL}_{LENGTH_DISTRIBUTION}_{ALIGNMENT_PROGRAM}_{new_file_index}.fasta"
    with open(true_msas_path/file_id, 'w') as f:
        f.write(sim_fasta_aligned)

    return sub_sim_msa_list, file_id

def compute_realigned_msa(unaligned_msa, file_id, realigned_msas_path):
    sim_fasta_unaligned = "".join([f'>{name}\n{seq}\n' for name, seq in zip(sparta_organism_order, unaligned_msa)])

    sim_fasta_unaligned = sim_fasta_unaligned.encode()
    with tempfile.NamedTemporaryFile(suffix='.fasta') as tempf:
        tempfile_path = pathlib.Path(tempf.name)
        tempf.write(sim_fasta_unaligned)
        tempf.seek(0)
        if ALIGNMENT_PROGRAM == "fastmafft":
            alignment_cline = MafftCommandline(input=tempf.name)
            realigned_msa, stderr = alignment_cline()
        if ALIGNMENT_PROGRAM == "mafft":
            alignment_cline = MafftCommandline(input=tempf.name, genafpair=True, maxiterate=1000)
            realigned_msa, stderr = alignment_cline()
        if ALIGNMENT_PROGRAM == "prank":
            alignment_cline = PrankCommandline(cmd="/groups/pupko/elyawygoda/bin/prank/prank/bin/prank",
                                            d=tempfile_path,
                                            o=tempfile_path)

            alignment_cline()
            with open(f"{tempfile_path}.best.fas" ,'r') as f:
                realigned_msa = f.read()
            os.unlink(f"{tempfile_path}.best.fas")

    with open(realigned_msas_path/file_id, 'w') as f:
        f.write(realigned_msa)


    realigned_msa = [s[s.index("\n"):].replace("\n","") for s in realigned_msa.split(">")[1:]]

    return realigned_msa

columns_sims = PARAMS_LIST + list(SUMSTATS_DEFINITION.values())
columns_realigned = [ f"RA_{col}" for col in SUMSTATS_DEFINITION.values()]
columns_realigned +=  ["INDEL_MODEL","LENGTH_DIST", "ALIGNMENT_PROGRAM", "ID"]
full_data_path = pathlib.Path(full_correction_path,"full_data.csv")

true_msas_path = pathlib.Path(full_correction_path, "simulated_msas")
try:
    print("creating true msas dir")
    os.mkdir(true_msas_path)
except:
    print("true msas folder exists already")

realigned_msas_path = pathlib.Path(full_correction_path, "realigned_msas")
try:
    print("creating realigned msas dir")
    os.mkdir(realigned_msas_path)
except:
    print("realigned folder exists already")

for msa, sim_stats in zip(msas_list, sum_stats_list):

    sim_msa, sim_file_id = add_substitutions(msa, true_msas_path)

    unaligned_msa = [seq.replace("-","") for seq in sim_msa]
    realigned_msa = compute_realigned_msa(unaligned_msa, sim_file_id, realigned_msas_path)

    file_index = sim_file_id.split("_")[-1].split(".")[0]
    realigned_stats = Sparta.Msa(realigned_msa).get_sum_stats()
    additional_info =  [INDEL_MODEL, LENGTH_DISTRIBUTION, ALIGNMENT_PROGRAM ,file_index]
    realigned_stats += additional_info

    sim_stats_df = pd.DataFrame([sim_stats], columns=columns_sims)
    realigned_stats_df = pd.DataFrame([realigned_stats], columns=columns_realigned)

    merged_df = pd.concat([sim_stats_df, realigned_stats_df], axis=1)
    if full_data_path.exists():
        merged_df.to_csv(full_data_path,mode="a", index=False, header=False)
    else:
        merged_df.to_csv(full_data_path, index=False)





