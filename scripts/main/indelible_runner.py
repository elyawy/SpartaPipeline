import os
import tempfile
import subprocess
from collections import OrderedDict

def get_indelible_config():
	indelible_config = OrderedDict()

	indelible_config["[TYPE]"] = 'AMINOACID 2'
	indelible_config["[MODEL]"] = 'modelname'
	indelible_config["[submodel]"] = 'WAG'
	indelible_config["[indelmodel]"] = 'POW 1.7 500'
	indelible_config["[indelrate]"] = '0.0'
	indelible_config["[rates]"] = ' 0.25 0.50 10'
	indelible_config["[statefreq]"] = ' 0.25  0.25  0.25  0.25'
	indelible_config["[TREE]"] = 'treename (A:0.1,B:0.1);'
	indelible_config["[PARTITIONS]"] = 'partitionname\n[treename modelname 3000]'
	indelible_config["[EVOLVE]"] = "partitionname 100 outputname" + " " # note: indlible requires space at last command.

	return indelible_config


def prepare_indelible_control_file(res_path, model_parameters):
	"""
	prepare indelible control file for simulating substitutions
	"""
	indelible_config = get_indelible_config()

	indelible_config["[TREE]"] = f'treename {model_parameters["tree"]}'
	indelible_config["[PARTITIONS]"] = f'partitionname\n[treename modelname {model_parameters["length"]}]'
	# note: do not change 'outputname1'
	indelible_config["[EVOLVE]"] = f'partitionname 1 outputname1' + " " # note: indlible requires space at last command.

	if model_parameters["mode"] == "amino":
		indelible_config["[TYPE]"] = 'AMINOACID 2'
		indelible_config["[submodel]"] = 'WAG'
		del indelible_config['[statefreq]']
		del indelible_config['[rates]']

	if model_parameters["mode"] == "nuc":
		indelible_config["[TYPE]"] = 'NUCLEOTIDE 2'
		
		if model_parameters["submodel"] == "GTR":
			gtr_params = ' '.join([f"{model_parameters['rates'][ind]:.9f}" for ind in range(5)])
			frequencies = ' '.join([f"{model_parameters['freq'][ind]:.6f}" for ind in range(4)])
			rates = f"{model_parameters['inv_prop']} {model_parameters['gamma_shape']} {model_parameters['gamma_cats']}"

			indelible_config["[submodel]"] = f'{model_parameters["submodel"]} {gtr_params}'
			indelible_config["[rates]"] = rates
			indelible_config["[statefreq]"] = frequencies

		if model_parameters["submodel"] == "JC":
			indelible_config["[submodel]"] = f'{model_parameters["submodel"]}'
			
			del indelible_config['[statefreq]']
			del indelible_config['[rates]']
	with open(os.path.join(res_path,'control.txt'),'w') as fout:
		for key in indelible_config:
			to_write = f'{key} {indelible_config[key]}\n'
			fout.write(to_write)


def IndelibleCommandline(model_params):
	"""
	runs indelible.
	Requires control.txt at res_path and indelible command
	"""

	origin_dir = os.getcwd()
	with tempfile.TemporaryDirectory() as tmpdirname:
		prepare_indelible_control_file(tmpdirname, model_params)
		os.chdir(tmpdirname)

		subprocess.run("indelible", stdout=subprocess.DEVNULL)
		indelible_msa_list = parse_indelible_output(tmpdirname)

		os.chdir(origin_dir)

	return indelible_msa_list

def parse_indelible_output(res_path):
	"""
	reads the output of indelible and parse it to list of msas
	"""
	with open(os.path.join(res_path,'outputname1.fas'),'r') as f:
		indelible_subs = f.read()

	indelible_msa_list = [s[s.index("\n"):].replace("\n","") for s in indelible_subs.split(">")[1:]]
	
	return indelible_msa_list



