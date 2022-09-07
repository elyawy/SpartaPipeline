from pathlib import Path
import argparse

from itertools import product

from elyawy.constants import length_distributions, indel_models

gr_parser = argparse.ArgumentParser(allow_abbrev=False)
gr_parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
gr_parser.add_argument('-db','--database', action='store',metavar="database to check", type=str, required=True)
args = gr_parser.parse_args()

MAIN_PATH = args.input
SELECTED_DB = args.database


all_models = product(length_distributions, indel_models)
all_models = ["_".join(model) for model in all_models]


results_path = Path(MAIN_PATH).resolve()


def check_results_folder(path_to_folder):
    results_status = []
    for model in all_models:
        current_path = path_to_folder / (f"full_data_{model}.pkl")
        current_correction_path = path_to_folder / f"correction/regressors_{model}"

        aggregate_check = str(int(current_path.exists())) + str(int(current_correction_path.exists()))
        results_status.append(aggregate_check)

    return '\t'.join(results_status)


test_summary = []
current_path = results_path / SELECTED_DB
for results in current_path.iterdir():
    test_summary.append(f"{results.stem}\t{check_results_folder(results)}")

print("data_name\t" + "\t".join(all_models))
test_summary = "\n".join(test_summary)
print(test_summary)