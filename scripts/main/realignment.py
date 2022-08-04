import tempfile, argparse, pathlib
from Bio.Align.Applications import MafftCommandline

gr_parser = argparse.ArgumentParser(allow_abbrev=False)
gr_parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)


args = gr_parser.parse_args()

MAIN_PATH = pathlib.Path(args.input)

MSA_PATH = list(filter(lambda x: x.suffix == ".fasta" ,MAIN_PATH.iterdir()))[0]



with open(MSA_PATH, 'r') as msaf:
    msa_string = msaf.read()

unaligned_fasta = msa_string.replace("-","")
unaligned_fasta = unaligned_fasta.encode()

with tempfile.NamedTemporaryFile(suffix='.fasta') as tempf:
    tempf.write(unaligned_fasta)
    tempf.seek(0)
    mafft_cline = MafftCommandline(input=tempf.name, genafpair=True, maxiterate=1000)
    realigned_msa, stderr = mafft_cline()


with open(pathlib.Path(MAIN_PATH,"realigned_msa.fasta"),'w') as f:
    f.write(realigned_msa)
