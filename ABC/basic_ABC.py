import sys
import numpy as np
import os
import config_reader



global one_out_m0
global one_out_m1
global one_out_a_param
one_out_m0 = None
one_out_m1 = None
one_out_a_param = [1, 2, 3, 19]

def extract_empriric(empiric_file, tree_file):
    """ receives a file of the empiric data and a tree file
    returns the summary statistics of the data"""
    return []

def simulation_generator(config, tree):
    """
    create simulation object
    :param config: dictionary of the configuration (priors etc)
    :param tree: tree obect
    :return: yields params and summary stats of generated simulation
    """
    yield [], []


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
    return params_set[1][2:4]


def extract_a_param(params_set):
    return [params_set[1][0]]


def extract_m1_params(params_set):
    return [len(params_set[0])] + params_set[1][4:6]


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
    sumstats = list(next(sumstat_generator)[1] for i in range(sim_number))  # does this works?
    ss_matrix = np.array(sumstats)
    ss_std = ss_matrix.std(axis=0)
    ss_std[ss_std == 0] = 1
    return 1 / ss_std


def filter_simulations_both_models_and_a(run_param_dict, ss_generator, weights):
    retained_sim = run_param_dict["retained_sim"]

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
    for params, full_ss in ss_generator:
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


def run_abc(emp_ss, conf, tree_file, out_path='None'):

    # create simulation generator
    tree = tree_file  # should either process to tree object or change generator to accept tree file
    sim_gen = simulation_generator()
    # calculate weights
    weights_all = compute_weights(sim_gen, conf["weight_size"])
    weights = {'m0': get_m0_ss(weights_all), 'a': get_a_param_ss(weights_all), 'm1': get_m1_ss(weights_all)}
    # filter simulations
    run_dict = {"retained_sim": conf["retained_sim"], "real_ss": emp_ss}
    m_0_param_list, a_param_list, m_1_param_list = filter_simulations_both_models_and_a(run_dict, sim_gen, weights)
    retained_dict = {'a': a_param_list, 'inversion': [m_0_param_list[i][0] for i in range(conf["retained_sim"])],
                     'translocation': [m_0_param_list[i][1] for i in range(conf["retained_sim"])],
                     'chr num': [m_1_param_list[i][0] for i in range(conf["retained_sim"])],
                     'fusion': [m_1_param_list[i][1] for i in range(conf["retained_sim"])],
                     'fission': [m_1_param_list[i][2] for i in range(conf["retained_sim"])]}
    if out_path is not 'None':
        for key in retained_dict:
            file_path = os.path.join(out_path, f"{key}_retained.txt")
            text = '\n'.join([str(x) for x in retained_dict[key]])
            with open(file_path, 'w') as f:
                f.write(text)
    return retained_dict



if __name__ == "__main__":
    empiric_file = sys.argv[1]
    tree_f = sys.argv[2]
    config_file = sys.argv[3]
    out_dir = sys.argv[4]  # can be None if we don't want to save the retained to file
    empiric_ss = extract_empriric(empiric_file, tree_f)
    prior_conf = config_reader.parse_conf(config_file)
    retained_dic = run_abc(empiric_ss, prior_conf, tree_f, out_dir)