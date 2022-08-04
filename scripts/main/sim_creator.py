import numpy as np
import config_reader

length_distribution_priors = {
    "zipf": {
        "insertion": [1.00000001, 1.93348850],
        "deletion": [1.00000001, 1.93348850]
    },
    "geometric": {
        "insertion": [0.02684840, 0.33333332],
        "deletion": [0.02684840, 0.33333332]
    },
    "poisson": {
        "insertion": [2.82143937, 20.0000001],
        "deletion": [2.82143937, 20.0000001]
    }
}

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
        
        insertion_lengths = [[i] for i in np.random.uniform(*self.len_prior_dict["insertion"], num_sims)]
        insertion_rates = np.random.uniform(*self.rate_prior_dict["insertion"], num_sims)

        if self.indel_model == "rim":
            deletion_lengths = [[i] for i in np.random.uniform(*self.len_prior_dict["deletion"], num_sims)]
            deletion_rates = np.random.uniform(*self.rate_prior_dict["deletion"], num_sims)
        elif self.indel_model == "sim":
            deletion_lengths = insertion_lengths
            deletion_rates = insertion_rates

        seq_length = np.random.randint(*self.sequence_length_prior, num_sims)

        len_dists = np.repeat(self.length_distribution, num_sims)

        return np.array([
                seq_length,
                len_dists,
                insertion_lengths, deletion_lengths,
                insertion_rates, deletion_rates,
               ], dtype=object).T
