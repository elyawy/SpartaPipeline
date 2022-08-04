import Sparta as _sp


class Msa:
    def __init__(self, msa):
        """
        Construct Msa object.

        Parameters
        ----------
            msa : str, Path, list[str]
                Path of the msa file in fasta format or list of aligned strings.
        """
        if isinstance(msa, (str)):
            self._msa = _sp.Msa(str(msa))
            self._msa.calc_stats()
        elif isinstance(msa, (list)):
            self._msa = _sp.Msa(msa)
            self._msa.calc_stats()
        elif isinstance(msa, (_sp.Msa)):
            self._msa = msa
            self._msa.calc_stats()

    def get_sum_stats(self):
        '''
        Return the summary statistics of the MSA, initialized to 0 except MSA length.
        '''
        return self._msa.get_sum_stats()
    

    def get_longest_seq_length(self):
        return self._msa.get_longest_seq_length()

    def get_shortest_seq_length(self):
        return self._msa.get_shortest_seq_length()

    def get_seq_with_indels(self):
        return self._msa.get_seq_with_indels()
        
    
    def get_seq_without_indels(self):
        return self._msa.get_seq_without_indels()

    def get_num_seq(self):
        return self._msa.get_num_seq()

    def print_msa(self):
        self._msa.print_msa()


class Simulator:
    def __init__(self, tree):
        """
        Construct Simulator object.
        By default the simulator is initialized with the following spec:
        root length: 100
        insertion rate: 0.03
        deletion rate: 0.03
        length distribution: zipf
        length parameter insertion: 1.5
        length parameter deletion: 1.5

        Parameters
        ----------
            tree : str, Path
                Path of the tree file in newick format
        """
        self.root_length = 100
        self.length_distribution = "zipf"
        self.length_parameter_insertion = [1.5]
        self.length_parameter_deletion = [1.5]
        self.insertion_rate = 0.03
        self.deletion_rate = 0.03

        self._sim = _sp.Sim(str(tree))
        self._sim.init_sim(self.root_length, self.length_distribution,
                          self.length_parameter_insertion, self.length_parameter_deletion,
                          self.insertion_rate, self.deletion_rate)


    def init_sim(self, root_length=100, insertion_rate=0.03, deletion_rate=0.03, length_distribution="zipf",
                 length_parameter_insertion=1.5, length_parameter_deletion=1.5):
        """
        Initalize the simulator with user given parameters.

        Parameters
        ----------
            root_length: int
            insertion_rate: float
            deletion_rate: float
            length_distribution: str
            length_parameter_insertion: float
            length_parameter_deletion: float
        """
        self.root_length = root_length
        self.length_distribution = length_distribution
        self.length_parameter_insertion = [length_parameter_insertion]
        self.length_parameter_deletion = [length_parameter_deletion]
        self.insertion_rate = insertion_rate
        self.deletion_rate = deletion_rate

        self._sim.init_sim(self.root_length, self.length_distribution,
                    self.length_parameter_insertion, self.length_parameter_deletion,
                    self.insertion_rate, self.deletion_rate)

    def set_root_length(self, root_length):
        """
        Set the root sequence length.

        Parameters
        ----------
            root_length: int
        """
        self.root_length = root_length

        self._sim.init_sim(self.root_length, self.length_distribution,
                    self.length_parameter_insertion, self.length_parameter_deletion,
                    self.insertion_rate, self.deletion_rate)

    def set_insertion_rate(self, insertion_rate):
        """
        Set the insertion rate.

        Parameters
        ----------
            insertion_rate: float
        """
        self.insertion_rate = insertion_rate

        self._sim.init_sim(self.root_length, self.length_distribution,
                    self.length_parameter_insertion, self.length_parameter_deletion,
                    self.insertion_rate, self.deletion_rate)
    
    def set_insertion_rate(self, deletion_rate):
        """
        Set the deletion rate.

        Parameters
        ----------
            deletion_rate: float
        """
        self.deletion_rate = deletion_rate

        self._sim.init_sim(self.root_length, self.length_distribution,
                    self.length_parameter_insertion, self.length_parameter_deletion,
                    self.insertion_rate, self.deletion_rate)
    
    def set_length_distribution(self, length_distribution):
        """
        Set the the length distribution used by the simulator.
        Options are: zipf, geometric, poisson.

        Parameters
        ----------
            length_distribution: str
        """
        self.length_distribution = length_distribution

        self._sim.init_sim(self.root_length, self.length_distribution,
                    self.length_parameter_insertion, self.length_parameter_deletion,
                    self.insertion_rate, self.deletion_rate)

    def set_length_parameter_insertion(self, length_parameter_insertion):
        """
        Set the insertion length parameter of the length distribution.

        Parameters
        ----------
            length_parameter_insertion: float
        """
        self.length_parameter_insertion = [length_parameter_insertion]

        self._sim.init_sim(self.root_length, self.length_distribution,
                    self.length_parameter_insertion, self.length_parameter_deletion,
                    self.insertion_rate, self.deletion_rate)

    def set_length_parameter_deletion(self, length_parameter_deletion):
        """
        Set the deletion length parameter of the length distribution.

        Parameters
        ----------
            length_parameter_deletion: float
        """
        self.length_parameter_deletion = [length_parameter_deletion]

        self._sim.init_sim(self.root_length, self.length_distribution,
                    self.length_parameter_insertion, self.length_parameter_deletion,
                    self.insertion_rate, self.deletion_rate)


    def set_seed(self, seed):
        """
        Set the random number generator seed.

        Parameters
        ----------
            seed: int
        """
        self._sim.set_seed(seed)

    def set_max_indel_size(self, max_size):
        """
        Set the truncation of the length distribution.

        Parameters
        ----------
            max_size: int
        """
        self._sim.set_max_indel_size(max_size)


    def __call__(self) -> Msa:
        '''
        Return
        ----------
            Msa: Simulated msa object.
        '''

        msa = self._sim.run_sim()
        return Msa(msa)
