import itertools
import random
import warnings

## takes output of has_mutations_by_genome_str
## and merges p/a (1/0) in genomes to individual level on a per-mutation basis
## (e.g. 101|001 -> 101, 010|001 -> 011)
## number of characters separated by '|' must be identical (function will not check)
def merge_genotype_str(s):
    phased = s.split('|')
    get_in_pos = lambda i: [genome[i] for genome in phased]
    output = ''.join('1' if '1' in get_in_pos(i) else '0' for i in range(len(phased[0])))
    return output

class SlimParseWarning(Warning):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

class ClassWithSubpopulation():
    def __init__(self, source, subpopulation_id):
        self.source = source
        self.subpopulation_id = subpopulation_id
    @property
    def subpopulation(self):
        return self.source.get_population(self.subpopulation_id)

class Population():
    def __init__(self, source, string = None, name = None, size = None, sex = None, sex_ratio = None):
        self.source = source
        self.name = name
        self.size = size
        self.sex = sex
        self.sex_ratio = sex_ratio
        self._individuals = {}
        self._genomes = {}
        if string is not None:
            self.parse_from_string(string)
        self.coerce_types()
    def coerce_types(self):
        ## make sure to coerce types
        if self.size is not None: self.size = int(self.size)
        if self.sex_ratio is not None: self.sex_ratio = float(self.sex_ratio)
        return
    def parse_from_string(self, string):
        ## string should look something like:
        ##   p1 50 H      (if hermaphrodites)
        ##   p2 50 S 0.5  (if separate sexes)
        split_str = string.split(' ')
        self.name, self.size, self.sex = split_str[:3]
        ## if population has separate sexes
        if self.sex == 'S':
            self.sex_ratio = split_str[3]
        self.coerce_types()
        return
    @property
    def individuals(self):
        return list(self._individuals.values())
    @property
    def individual_ids(self):
        return list(self._individuals.keys())
    def add_individual(self, individual):
        self._individuals[individual.id] = individual
        return
    def get_individual(self, id):
        return self._individuals.get(id, None)
    def get_individuals(self, *ids):
        return [self.get_individual(id) for id in ids]
    @property
    def genomes(self):
        return list(self._genomes.values())
    def add_genome(self, genome):
        self._genomes[genome.id] = genome
        return
    def get_genome(self, id):
        return self._genomes.get(id, None)
    def get_genomes(self, *ids):
        return [self.get_genome(id) for id in ids]
    @property
    def mutations(self):
        return list(set(itertools.chain(*[genome.mutations for genome in self.genomes])))
    @property
    def mutation_ids(self):
        return list(set(itertools.chain(*[genome.mutation_ids for genome in self.genomes])))
    def sample_individuals(self, n = 1, condition = lambda x:True, replacement = False):
        if replacement:
            ## sample with replacement
            ## check if there are individuals that satisfy condition(s)
            individuals_satisfying_condition = [indv for indv in self.individuals if condition(indv)]
            if len(individuals_satisfying_condition) == 0:
                ## if there are no individuals that pass the condition(s)
                warnings.warn("No individuals satisfy the desired condition(s). Returning empty list.",
                              SlimParseWarning)
                return []
            else:
                ## return list of individuals
                return random.choices(individuals_satisfying_condition, k = n)
        else:
            ## sample without replacement
            individuals = self.individuals
            n_individuals = len(individuals)
            ## randomly order individuals
            sample_order = random.sample(individuals, n_individuals)
            ## reduce n if necessary (i.e. n > num individuals)
            if n > n_individuals:
                n = n_individuals
                warnings.warn(("n cannot be larger than the number of individuals if replacement=False."
                               f" Setting n to number of individuals in simulation ({n_individuals})."),
                              SlimParseWarning)
            ## select first n individuals that satisfy condition(s)
            sample_individuals = []
            i = 0
            while len(sample_individuals) < n and i < n_individuals:
                individual = sample_order[i]
                if condition(individual):
                    sample_individuals.append(individual)
                i += 1
            ## warn if number of individuals that satisfy condition(s) < n
            if len(sample_individuals) < n:
                warnings.warn((f"Warning: There are fewer individuals with desired condition(s) than n ({n})."
                               f" Returning {len(sample_individuals)} sampled individuals."),
                              SlimParseWarning)
            return sample_individuals

class Mutation(ClassWithSubpopulation):
    def __init__(self, source, string = None, id = None, id_run = None, type = None,
                 position = None, selection_coefficient = None, dominance_coefficient = None,
                 origin_subpopulation_id = None, origin_tick = None, prevalence = None, base = None):
        super().__init__(source, origin_subpopulation_id)
        self.id = id
        self.id_run = id_run
        self.type = type
        self.position = position
        self.selection_coefficient = selection_coefficient
        self.dominance_coefficient = dominance_coefficient
        self.origin_tick = origin_tick
        self.prevalance = prevalence
        self.base = base
        if string is not None:
            self.parse_from_string(string)
        self.coerce_types()
    @property
    def origin_subpopulation(self):
        return self.subpopulation
    @property
    def origin_subpopulation_id(self):
        return self.subpopulation_id
    @origin_subpopulation_id.setter
    def origin_subpopulation_id(self, val):
        self.subpopulation_id = val
        return
    def coerce_types(self):
        ## make sure to coerce types
        if self.id is not None: self.id = int(self.id)
        if self.id_run is not None: self.id_run = int(self.id_run)
        if self.position is not None: self.position = int(self.position)
        if self.selection_coefficient is not None: self.selection_coefficient = float(self.selection_coefficient)
        if self.dominance_coefficient is not None: self.dominance_coefficient = float(self.dominance_coefficient)
        if self.origin_tick is not None: self.origin_tick = int(self.origin_tick)
        if self.prevalence is not None: self.prevalence = float(self.prevalence)
        return
    def parse_from_string(self, string):
        ## string should look something like:
        ##   10 387752 m1 1308 0 0.5 p2 9404 130    (if not nucleotide-based)
        ##   47 387966 m2 5994 0 0.5 p1 9415 130 A  (if nucleotide-based)
        split_str = string.split(' ')
        self.id, self.id_run, self.type, self.position = split_str[:4]
        self.selection_coefficient, self.dominance_coefficient = split_str[4:6]
        self.origin_subpopulation_id, self.origin_tick, self.prevalence = split_str[6:9]
        if len(split_str) > 9:
            self.base = split_str[9]
        self.coerce_types()
        return
    def is_nucleotide_based(self):
        return (self.base is not None)
    def is_mutation_type(self, mutation_type: str):
        return self.type == mutation_type

class Individual(ClassWithSubpopulation):
    def __init__(self, source, string = None, id = None, pedigree_id = None, sex = None,
                 genome_1_id = None, genome_2_id = None, position = None, subpopulation_id = None):
        self.id = id
        subpopulation_id = (subpopulation_id if subpopulation_id is not None else
                            None if self.id is None
                            else self.id.split(':')[0]) ## extract subpopulation from id if id is not None
        super().__init__(source, subpopulation_id)
        self.pedigree_id = pedigree_id ## present only if pedigreeIDs=T (specified in OutputFull obj)
        self.sex = sex
        self.genome_ids = [genome_1_id, genome_2_id] ## list of genome ids (str)
        self.position = position ## present only if spatialPositions=T (specified in OutputFull obj)
        if string is not None:
            self.parse_from_string(string)
        self.coerce_types()
    def coerce_types(self):
        if self.position is not None:
            self.position = [float(pos) for pos in self.position]
        return
    def parse_from_string(self, string):
        split_str = string.split(' ')
        self.id = split_str[0]
        self.subpopulation_id = self.id.split(':')[0]
        if self.source.pedigreeIDs:
            self.pedigree_id = split_str[1]
            i = 2
        else:
            i = 1
        self.sex = split_str[i]
        self.genome_ids = split_str[i+1:i+3]
        if self.source.spatialPositions:
            self.position = split_str[i+3:]
        else:
            self.position = None
        self.coerce_types()
        return
    @property
    def genomes(self):
        ## retrieves Genome objects based on genome ids
        return self.subpopulation.get_genomes(*self.genome_ids)
    @property
    def x(self):
        return None if self.position is None else self.position[0]
    @property
    def y(self):
        return None if self.position is None else self.position[1]
    @property
    def z(self):
        return None if (self.position is None or len(self.position) < 3) else self.position[2]
    @property
    def mutations(self):
        return list(set(itertools.chain(*[genome.mutations for genome in self.genomes])))
    @property
    def mutation_ids(self):
        return list(set(itertools.chain(*[genome.mutation_ids for genome in self.genomes])))
    @property
    def mutation_types(self):
        return list(set(itertools.chain(*[genome.mutation_types for genome in self.genomes])))
    ## returns list of T/F, where each element is presence of a mutation/mutation type
    ## in a given genome in this individual
    def mutation_in_genomes(self, mutation):
        ## accepts Mutation obj, int of mutation ID, and str of mutation type
        return [genome.has_mutation(mutation) for genome in self.genomes]
    # ## returns list of T/F, where each element is presence of all mutation/mutation types in 'mutations' inpt
    # ## in a given genome in this individual
    # def all_mutations_in_genomes(self, *mutations):
    #     ## accepts Mutation obj, int of mutation ID, and str of mutation type
    #     return [genome.has_all_mutations(*mutations) for genome in self.genomes]
    # def any_mutations_in_genomes(self, *mutations):
    #     ## accepts Mutation obj, int of mutation ID, and str of mutation type
    #     return [genome.has_any_mutations(*mutations) for genome in self.genomes]
    def has_mutation(self, mutation):
        ## accepts Mutation obj, int of mutation ID, and str of mutation type
        return any(self.mutation_in_genomes(mutation))
    def is_homozygous(self, mutation):
        ## accepts Mutation obj, int of mutation ID, and str of mutation type
        return all(self.mutation_in_genomes(mutation))
    def is_heterozygous(self, mutation):
        ## accepts Mutation obj, int of mutation ID, and str of mutation type
        return sum(self.mutation_in_genomes(mutation)) == 1
    def has_all_mutations(self, *mutations):
        return all(self.has_mutation(mutation) for mutation in mutations)
    def has_any_mutations(self, *mutations):
        return any(self.has_mutation(mutation) for mutation in mutations)
    ## returns list of T/F, where each element is presence of all mutation/mutation types in 'mutations' inpt
    ## in this individual
    def has_mutations(self, *mutations):
        return [self.has_mutation(mutation) for mutation in mutations]
    ## returns [[T/F of mutations P/A in genome 1], [T/F of mutations P/A in genome 2]]
    def has_mutations_by_genome(self, *mutations, sort = False):
        output = [genome.has_mutations(*mutations) for genome in self.genomes]
        if sort: return sorted(output)
        else: return output
    ## converts list of T/F from self.has_mutations to str of 1/0 (where 1=T)
    def has_mutations_str(self, *mutations):
        return ''.join('1' if present else '0' for present in self.has_mutations(*mutations))
    ## returns <str of mutation P/A in 1/0 of genome 1>|<str of mutation P/A in 1/0 of genome 2>
    def has_mutations_by_genome_str(self, *mutations, sort = False):
        output = [genome.has_mutations_str(*mutations) for genome in self.genomes]
        if sort: output.sort()
        return '|'.join(output)

class Genome(ClassWithSubpopulation):
    def __init__(self, source, string = None, id = None, type = None, mutation_ids = None,
                 subpopulation_id = None):
        self.id = id
        subpopulation_id = (subpopulation_id if subpopulation_id is not None else
                            None if self.id is None
                            else self.id.split(':')[0]) ## extract subpopulation from id if id is not None
        super().__init__(source, subpopulation_id)
        self.type = type
        self._mutation_ids_set = None
        self._mutation_ids = None
        self.mutation_ids = mutation_ids
        if string is not None:
            self.parse_from_string(string)
        self.coerce_types()
    @property
    def mutation_ids(self):
        return self._mutation_ids
    @mutation_ids.setter
    def mutation_ids(self, new_ids):
        if new_ids is not None:
            self._mutation_ids = [int(id) for id in new_ids]
            self._mutation_ids_set = set(self._mutation_ids)
        else:
            self._mutation_ids = new_ids
            self._mutation_ids_set = None
        return
    def coerce_types(self):
        if self.mutation_ids is not None:
            self.mutation_ids = [int(id) for id in self.mutation_ids]
        return
    def parse_from_string(self, string):
        split_str = string.split(' ')
        self.id, self.type = split_str[:2]
        self.subpopulation_id = self.id.split(':')[0]
        if len(split_str) > 2 and split_str[2] == "<null>":
            self.mutation_ids = []
        else:
            self.mutation_ids = split_str[2:]
        self.coerce_types()
        return
    @property
    def mutations(self):
        ## list of Mutation objects in this genome
        return [self.source.get_mutation(id) for id in self.mutation_ids]
    @property
    def mutation_types(self):
        ## list of mutation types (str) in this genome (after applying set)
        return list(set([mutation.type for mutation in self.mutations]))
    ## returns true if at least one instance of mutation/mutation type is in this genome
    def has_mutation(self, mutation):
        ## accepts Mutation obj, int of mutation ID, and str of mutation type
        if isinstance(mutation, str): ## if input is str of mutation type
            return self.has_mutation_type(mutation)
        else: ## elif input is int of mutation ID or Mutation obj
            mutation_id = mutation if isinstance(mutation, int) else mutation.id
            return mutation_id in self._mutation_ids_set
    ## returns true if at least one mutation of mutation_type is in this genome
    def has_mutation_type(self, mutation_type):
        return any(mutation.is_mutation_type(mutation_type) for mutation in self.mutations)
    ## returns true if all mutations are in this genome
    def has_all_mutations(self, *mutations):
        ## accepts mix of Mutation obj, int of mutation ID, and str of mutation type
        return all(self.has_mutation(mutation) for mutation in mutations)
    ## returns true if one or more of the mutations are in this genome
    def has_any_mutations(self, *mutations):
        ## accepts mix of Mutation obj, int of mutation ID, and str of mutation type
        return any(self.has_mutation(mutation) for mutation in mutations)
    ## returns list of T/F, where each element is presence of all mutation/mutation types in 'mutations' inpt
    ## in this genome
    def has_mutations(self, *mutations):
        return [self.has_mutation(mutation) for mutation in mutations]
    ## converts list of T/F from self.has_mutations to str of 1/0 (where 1=T)
    def has_mutations_str(self, *mutations):
        return ''.join('1' if present else '0' for present in self.has_mutations(*mutations))

class OutputFull():
    def __init__(self, fname = None, string = None, pedigreeIDs = False, spatialPositions = False):
        self.fname = None
        self.string = None
        self.version = None
        self.pedigreeIDs = pedigreeIDs
        self.spatialPositions = spatialPositions
        self.tick = None
        self.cycle = None
        self._populations = {}
        self._mutations = {}
        self.parse(fname = fname, string = string)
    def parse(self, fname = None, string = None):
        if fname is not None:
            self.parse_from_file(fname)
        elif string is not None:
            self.parse_from_string(string)
        else:
            warnings.warn(f"fname or string is required for {self.__class__.__name__}.parse()",
                          SlimParseWarning)
        return
    def parse_from_file(self, fname):
        with open(fname, 'r') as f:
            string = f.read()
            self.parse_from_string(string)
        return
    def parse_from_string(self, string):
        lines = string.splitlines()
        n_lines = len(lines)
        ## get some metadata
        tick, cycle = lines[0].split(' ')[1:3]
        self.tick = int(tick)
        self.cycle = int(cycle)
        if "Version:" in lines[1]:
            self.version = lines[1].replace("Version: ", '')
            i = 3 ## row index of first populations data
        else:
            i = 2 ## row index of first populations data
        ## function to iterate through lines
        def iterate(start, continuation_condition, action):
            i = start
            while continuation_condition(i):
                action(i)
                i += 1
            return i
        ## parse populations
        i = iterate(i, lambda i: lines[i] != "Mutations:", lambda i: self.parse_population(lines[i]))
        ## parse mutations
        i += 1
        i = iterate(i, lambda i: lines[i] != "Individuals:", lambda i: self.parse_mutation(lines[i]))
        any_nucleotide_based_muts = any(mut.is_nucleotide_based() for mut in self.mutations)
        ## parse individuals
        i += 1
        i = iterate(i, lambda i: lines[i] != "Genomes:", lambda i: self.parse_individual(lines[i]))
        ## parse genomes (continuation condition depends on whether 'Ancestral sequence:' block exists)
        if any_nucleotide_based_muts:
            is_genome_entry = lambda i: lines[i] != "Ancestral sequence:"
        else:
            is_genome_entry = lambda i: i < n_lines
        i += 1
        i = iterate(i, is_genome_entry, lambda i: self.parse_genome(lines[i]))
        ## parse ancestral sequences
        i += 1
        i = iterate(i, lambda i: i < n_lines, lambda i: self.parse_ancestral_sequence(lines[i]))
        return
    def add_population(self, population):
        self._populations[population.name] = population
    def add_mutation(self, mutation):
        self._mutations[mutation.id] = mutation
    def get_population(self, population_name):
        return self._populations[population_name]
    def get_mutation(self, mutation_id):
        return self._mutations[mutation_id]
    @property
    def populations(self):
        return list(self._populations.values())
    @property
    def mutations(self):
        return list(self._mutations.values())
    @property
    def population_names(self):
        return list(self._populations.keys())
    @property
    def mutation_ids(self):
        return list(self._mutations.keys())
    def parse_population(self, string):
        pop_name = string[:string.index(' ')]
        self.add_population(Population(self, string = string))
    def parse_mutation(self, string):
        mutation_id = string[:string.index(' ')]
        self.add_mutation(Mutation(self, string = string))
    def parse_individual(self, string):
        indv_id = string[:string.index(' ')]
        pop_name = indv_id.split(':')[0]
        self.get_population(pop_name).add_individual(Individual(self, string = string))
    def parse_genome(self, string):
        genome_id = string[:string.index(' ')]
        pop_name = genome_id.split(':')[0]
        self.get_population(pop_name).add_genome(Genome(self, string = string))
    @property
    def individuals(self):
        return [indv for population in self.populations for indv in population.individuals]
    @property
    def genomes(self):
        return [genome for population in self.populations for genome in population.genomes]
    def sample_mutations(self, n = 1, condition = lambda x:True, replacement = False):
        if replacement:
            ## sample with replacement
            ## check if there are mutations that satisfy condition(s)
            mutations_satisfying_condition = [mut for mut in self.mutations if condition(mut)]
            if len(mutations_satisfying_condition) == 0:
                ## if there are no mutations that pass the condition(s)
                warnings.warn("No mutations satisfy the desired condition(s). Returning empty list.",
                              SlimParseWarning)
                return []
            else:
                ## return list of mutations
                return random.choices(mutations_satisfying_condition, k = n)
        else:
            ## sample without replacement
            ## randomly order mutations
            mutations = self.mutations
            n_mutations = len(mutations)
            sample_order = random.sample(mutations, n_mutations)
            ## reduce n if necessary (i.e. n > num mutations)
            if n > n_mutations:
                warnings.warn(("n cannot be larger than the number of mutations if replacement=False."
                               f" Setting n to number of mutations in simulation ({n_mutations})."),
                              SlimParseWarning)
                n = n_mutations
            ## select first n mutations that satisfy condition(s)
            sample_mutations = []
            i = 0
            while len(sample_mutations) < n and i < n_mutations:
                mutation = sample_order[i]
                if condition(mutation):
                    sample_mutations.append(mutation)
                i += 1
            ## warn if number of mutations that satisfy condition(s) < n
            if len(sample_mutations) < n:
                warnings.warn((f"Warning: There are fewer mutations with desired condition(s) than n ({n})."
                               f" Returning {len(sample_mutations)} sampled mutations."),
                              SlimParseWarning)
            return sample_mutations


# ## test
# fname = "/mnt/chaelab/rachelle/scd/results/recipe_run_20230911_sample/scd_1_1/slim_out/scd_1_1752856844113_30.txt"
# sim = OutputFull(fname = fname)
# sim.genomes[2].mutations
# sim.sample_mutations(2)
# sim.get_population("p1").sample_individuals(2)

# ## get m2,m3,m4 P/A for all individuals (phased)
# import sys
# sys.path.append("/mnt/chaelab/rachelle/src")
# from basics import get_count_dict

# ## get p/a of entire population
# pa_m234 = [indv.has_mutations_by_genome_str("m2", "m3", "m4") for indv in sim.get_population("p1").individuals]
# pa_m234_unphased = [merge_genotype_str(pa) for pa in pa_m234]
# pa_m234_count = get_count_dict(pa_m234)
# pa_m234_unphased_count = get_count_dict(pa_m234_unphased)

# ## select neutral mutations (mutationType = "m1") to act as control
# p1_size = sim.get_population("p1").size
# min_freq, max_freq = 0.05*p1_size, 0.95*p1_size
# control_muts = sim.sample_mutations(5, condition = lambda mut: mut.type == "m1" and mut.prevalence >= min_freq and mut.prevalence <= max_freq)

# ## get p/a of subsample of population for m2-neutral mutations (this should be repeated for m2-m3 too)
# ## this section should be run at least 10 times per output file (and per n_indvs_to_sample (possibly 50, 100, 200?)) to get an idea of sampling size effect
# n_indvs_to_sample = 50
# sample_indvs = sim.get_population("p1").sample_individuals(n_indvs_to_sample)
# pa_m2x4_unphased = {mut.id: [indv.has_mutations_str("m2", mut.id, "m4")
#                              for indv in sample_indvs]
#                     for mut in control_muts}
# pa_m2x4_unphased_count = {mut_id: get_count_dict(data) for mut_id, data in pa_m2x4_unphased.items()}
