import argparse
import numpy as np
import os
import GenomeRearrangement as GR
import config_reader

gr_parser = argparse.ArgumentParser(allow_abbrev=False)
gr_parser.add_argument('-i','--input', action='store',metavar="Empirical genome", type=str, required=True)
gr_parser.add_argument('-c','--config', action='store',metavar="Simulation config" , type=str, required=True)
gr_parser.add_argument('-t','--tree', action='store', metavar="Newick tree", type=str, required=True)
gr_parser.add_argument('-o','--output', action='store', metavar="Output directory", type=str, required=False)



global one_out_m0
global one_out_m1
global one_out_a_param
one_out_m0 = None
one_out_m1 = None
one_out_a_param = [1, 2, 3, 19]

def extract_empriric(empiric_file, tree_file):
    """ receives a file of the empiric data and a tree file
    returns the summary statistics of the data"""
    sumstats = GR.genomes(empiric_file, tree_file).get_sum_stats()
    return sumstats


def params_generator(config):
    """
    generate params list
    :param config: dictionary of the configuration (priors etc)
    :return: yields parameters from given priors
    """
    def get_chromosomes(genome_size, numChromosomes):
        counts = np.random.multinomial(int(genome_size), [1.0/numChromosomes]*numChromosomes)
        return counts.tolist()
    def get_params():
        params = {
            "Chromosomes": get_chromosomes(config["GenomeSize"](), config["numChromosomes"]()),
            "AParam": config["AParam"](),
            "maxBlockSize": int(config["maxBlockSize"]()),
            "InvertRate": config["InvertRate"](),
            "TranslocateRatio": config["TranslocateRatio"](),
            "FusionRate": config["FusionRate"](),
            "FissionRate": config["FissionRate"](),
            "DuplicationRate": config["DuplicationRate"](),
            "LossRate":  config["LossRate"](),
            "randRootAparam":  config["randRootAparam"]()
        }
        return list(params.values())
    return get_params

def get_simulation_generator(config, tree_file):
    """
    create simulation object
    :param config: dictionary of the configuration (priors etc)
    :param tree: tree file
    :return: yields params and summary stats of generated simulation
    """
    params = params_generator(config)
    simulator = GR.Sim(tree_file)
    def simulation_generator(sim_num):
        for i in range(sim_num):
            sim_params = params()
            simulator.init_sim(*sim_params)
            simulator.set_seed(int(config['seed']()))

            genome_object = GR.genomes(simulator.run_sim(), simulator.get_tree())

            yield sim_params, genome_object.get_sum_stats()
    return simulation_generator

def get_m0_ss(full_ss):
    return leave_out(np.array(full_ss[0] + full_ss[1]), one_out_m0)


def get_a_param_ss(full_ss):
    return leave_out(np.array(full_ss[0] + full_ss[1]), one_out_a_param)


def get_m1_ss(full_ss):
    return leave_out(np.array(full_ss[2]), one_out_m1)


def leave_out(ss, index):
    # remove specified summary statistic from vector
    if index is None:
        return ss
    return np.delete(ss, index)


def extract_m0_params(params_set):
    return [params_set[1][0]] + params_set[1][2:4]


def extract_m0_params_no_a(params_set):
    return params_set[3:5]


def extract_a_param(params_set):
    return params_set[1]


def extract_m1_params(params_set):
    return [len(params_set[0])] + params_set[5:7]


def compute_distance(vec1, vec2, weight_vec):
    # print(f"v1 {vec1}")
    # print(f"v2 {vec2}")
    assert len(vec1) == len(vec2)
    # print(weight_vec)
    assert len(vec1) == len(weight_vec)
    return (weight_vec*(vec2 - vec1)**2).sum()  # squared distance


def compute_weights(sumstat_generator, sim_number):
    # receives a generator of summary stats
    # calculates the weights and returns the weight vector
    sumstats_list = [[],[],[]]
    for params, stats in sumstat_generator(sim_number):
        for j in range(3):
            sumstats_list[j].append(stats[j])
    
    for j in range(3):
        ss_matrix = np.array(sumstats_list[j])
        ss_std = ss_matrix.std(axis=0)
        ss_std[ss_std == 0] = 1
        sumstats_list[j] = (1/ss_std).tolist()
    return sumstats_list


def filter_simulations_both_models_and_a(run_param_dict, ss_generator, weights):
    retained_sim = run_param_dict["retained_sim_num"]

    def update_retained(ret_list, dist, params):
        if len(ret_list) < retained_sim or (dist, params) < max(ret_list):
            if len(ret_list) == retained_sim:
                ret_list.pop()
            ret_list.append((dist, params))
            ret_list.sort()
            return

    def prep_to_return(ret_list):
        # print(f" retained list size {len(ret_list)}")
        return [pair[1] for pair in ret_list]
    """ filter for both models and a param"""
    real_ss_m0 = get_m0_ss(run_param_dict["real_ss"])
    real_ss_m1 = get_m1_ss(run_param_dict["real_ss"])
    real_ss_a = get_a_param_ss(run_param_dict["real_ss"])
    dist_ind_list_m0 = []  # will hold the distance, index couples
    dist_ind_list_a = []  # will hold the distance, index couples
    dist_ind_list_m1 = []  # will hold the distance, index couples
    for params, full_ss in ss_generator(run_param_dict['sim_num']):
        sim_ss_m0 = get_m0_ss(full_ss)
        sim_ss_m1 = get_m1_ss(full_ss)
        sim_ss_a = get_a_param_ss(full_ss)
        distance_m0 = compute_distance(real_ss_m0, sim_ss_m0, weights["m0"])
        distance_a = compute_distance(real_ss_a, sim_ss_a, weights["a"])
        distance_m1 = compute_distance(real_ss_m1, sim_ss_m1, weights["m1"])
        update_retained(dist_ind_list_m0, distance_m0, extract_m0_params_no_a(params))
        update_retained(dist_ind_list_a, distance_a,  extract_a_param(params))
        update_retained(dist_ind_list_m1, distance_m1, extract_m1_params(params))
    return prep_to_return(dist_ind_list_m0), prep_to_return(dist_ind_list_a), prep_to_return(dist_ind_list_m1)


def run_abc(emp_ss, conf, tree_file, out_path=None):

    # create simulation generator
    sim_gen = get_simulation_generator(conf, tree_file=tree_file)
    # calculate weights
    weights_all = compute_weights(sim_gen, int(conf["weight_size"]()))
    weights = {'m0': get_m0_ss(weights_all), 'a': get_a_param_ss(weights_all), 'm1': get_m1_ss(weights_all)}
    # filter simulations
    run_dict = {"retained_sim_num": int(conf["retained_sim_num"]()), "sim_num": int(conf["sim_num"]()), "real_ss": emp_ss}
    m_0_param_list, a_param_list, m_1_param_list = filter_simulations_both_models_and_a(run_dict, sim_gen, weights)
    retained_dict = {'a': a_param_list, 'inversion': [m_0_param_list[i][0] for i in range(run_dict["retained_sim_num"])],
                     'translocation': [m_0_param_list[i][1] for i in range(run_dict["retained_sim_num"])],
                     'chr num': [m_1_param_list[i][0] for i in range(run_dict["retained_sim_num"])],
                     'fusion': [m_1_param_list[i][1] for i in range(run_dict["retained_sim_num"])],
                     'fission': [m_1_param_list[i][2] for i in range(run_dict["retained_sim_num"])]}
    if out_path is not None:
        for key in retained_dict:
            file_path = os.path.join(out_path, f"{key}_retained.txt")
            text = '\n'.join([str(x) for x in retained_dict[key]])
            with open(file_path, 'w') as f:
                f.write(text)
    return retained_dict



if __name__ == "__main__":
    args = gr_parser.parse_args()
    empiric_ss = extract_empriric(args.input, args.tree)
    prior_conf = config_reader.parse_conf(args.config)
    np.random.seed(int(prior_conf['seed']()))
    retained_dic = run_abc(empiric_ss, prior_conf, args.tree, args.output)