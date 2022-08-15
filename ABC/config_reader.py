import numpy as np

def const_func(prior):
    def get_const_number():
        return prior[0]
    return get_const_number

def uniform_int_func(prior):
    def get_uniform_int():
        return np.random.randint(*prior)
    return get_uniform_int

def uniform_float_func(prior):
    def get_uniform_float():
        return np.random.uniform(*prior)
    return get_uniform_float

def log_normal_func(prior):
    def get_log_normal():
        return np.random.lognormal(*prior)
    return get_log_normal

def gamma_func(prior):
    def get_gamma():
        return np.random.gamma(*prior)
    return get_gamma

def poisson_func(prior):
    def get_poisson():
        return np.random.poisson(*prior)
    return get_poisson

np_functions = {
    "uniform_int": uniform_int_func,
    "uniform_float": uniform_float_func,
    "log_normal": log_normal_func,
    "gamma": gamma_func,
    "poisson": poisson_func,
    "const": const_func
}

default_params = {
    "sim_num": lambda: 10,
    "seed": lambda: 1,
    "tree_path": lambda: ".",
    "output_dir_path": lambda: ".",
    "GenomeSize": lambda: 1000,
    "numChromosomes": (np.random.randint ,3,10),
    "AParam": (np.random.uniform, 1.01,1.99),
    "maxBlockSize": lambda: 50,
    "InvertRate": (np.random.uniform, 0.01,0.99),
    "TranslocateRatio": (np.random.uniform, 0.01,0.99),
    "FusionRate": (np.random.uniform, 0.01,0.99),
    "FissionRate": (np.random.uniform, 0.01,0.99),
    "DuplicationRate": lambda: 0.0,
    "LossRate": lambda: 0.0,
    "randRootAparam": lambda: 1.5,
    "retained_sim_num":lambda: 50
}



def parse_conf(path_to_conf_file):
    params = default_params.copy()
    with open(path_to_conf_file,'r') as f:
        for line in f:
            name, values = line.strip().split(":")
            values = values.strip().split(",")
            prior = [float(x.strip()) for x in values[1:]]
            params[name] = np_functions[values[0].strip()](prior)

    return params