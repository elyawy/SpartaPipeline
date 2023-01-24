import pickle, pathlib, os
import numpy as np

import config_reader
from elyawy.sparta import Simulator



current_path = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
priors_path = (current_path / 'getting_priors/priors.pickle').resolve()
length_distribution_priors = {}
if priors_path.exists():
    with open(priors_path, 'rb') as handle:
        length_distribution_priors = pickle.load(handle)
else:
    from getting_priors import get_means

    length_distribution_priors = {
        "zipf": {
            "insertion": sorted(get_means.final_priors["zipf"]),
            "deletion": sorted(get_means.final_priors["zipf"])
        },
        "geometric": {
            "insertion": sorted(get_means.final_priors["geometric"]),
            "deletion": sorted(get_means.final_priors["geometric"])
        },
        "poisson": {
            "insertion": sorted(get_means.final_priors["poisson"]),
            "deletion": sorted(get_means.final_priors["poisson"])
        }
    }

    with open(priors_path, 'wb') as handle:
        pickle.dump(length_distribution_priors, handle, protocol=pickle.HIGHEST_PROTOCOL)



class SimConfig:
    def __init__(self, conf_file=None,
                       len_dist="zipf",
                       rate_priors=[[0.01,0.05],[0.01,0.05]],
                       seq_lengths=[100,500],
                       indel_model="sim"):
        self.seed = 1
        if conf_file is None:
            self.indel_model = indel_model
            self.length_distribution = len_dist

            self.len_prior_dict = length_distribution_priors[len_dist]

            self.rate_prior_dict = {
                "insertion": rate_priors[0],
                "deletion": rate_priors[1],
            }
        else:
            configuration = config_reader.parse_conf(conf_file)
            self.seed = int(configuration["seed"])
            self.indel_model = configuration["indel_model"]
            self.length_distribution = configuration["length_distribution"]
            self.len_prior_dict = length_distribution_priors[self.length_distribution]
            self.rate_prior_dict = {
                "insertion": rate_priors[0],
                "deletion": rate_priors[1],
            }
        self.sequence_length_prior = [seq_lengths[0]*0.8, seq_lengths[1]*1.1]




    def get_random_sim(self, num_sims):
        
        insertion_lengths = [i for i in np.random.uniform(*self.len_prior_dict["insertion"], num_sims)]
        insertion_rates = np.random.uniform(*self.rate_prior_dict["insertion"], num_sims)

        if self.indel_model == "rim":
            deletion_lengths = [i for i in np.random.uniform(*self.len_prior_dict["deletion"], num_sims)]
            deletion_rates = np.random.uniform(*self.rate_prior_dict["deletion"], num_sims)
        elif self.indel_model == "sim":
            deletion_lengths = insertion_lengths
            deletion_rates = insertion_rates

        root_lengths = np.random.randint(*self.sequence_length_prior, num_sims)

        length_distribution = np.repeat(self.length_distribution, num_sims)

        #     root_lengths: int
        #     insertion_rates: float
        #     deletion_rates: float
        #     length_distribution: str
        #     insertion_lengths: float
        #     deletion_lengths: float

        return np.array([
                root_lengths,
                insertion_rates,
                deletion_rates,
                length_distribution,
                insertion_lengths, deletion_lengths
                ], dtype=object).T
