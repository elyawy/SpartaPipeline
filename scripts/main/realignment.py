import tempfile, argparse, pathlib
from aligner_interface import Aligner

_parser = argparse.ArgumentParser(allow_abbrev=False)
_parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
_parser.add_argument('-a','--aligner', action='store',metavar="Alignment program to use" , type=str, required=True)


args = _parser.parse_args()

MAIN_PATH = pathlib.Path(args.input)
ALIGNER = Aligner(args.aligner)

MSA_PATH = list(filter(lambda x: x.suffix == ".fasta" ,MAIN_PATH.iterdir()))[0]



with open(MSA_PATH, 'r') as msaf:
    msa_string = msaf.read()

unaligned_fasta = msa_string.replace("-","")
unaligned_fasta = unaligned_fasta.encode()

with tempfile.NamedTemporaryFile(suffix='.fasta') as tempf:
    tempf.write(unaligned_fasta)
    tempf.seek(0)
    ALIGNER.set_input_file(tempf.name)
    realigned_msa = ALIGNER.get_realigned_msa()


with open(pathlib.Path(MAIN_PATH,f"{args.aligner}_realigned_msa.fasta"),'w') as f:
    f.write(realigned_msa)
