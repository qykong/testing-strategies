import argparse
import sys
import os
import random
import json
from pyscripts.utils import *

DEFAULT_TRACING_PARAMS_FILE = 'parameters/default_tracing_parameters.json'
DEFAULT_POPULATION_PARAMS_FILE = 'parameters/default_population_parameters.json'

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

class NpEncoder(json.JSONEncoder):
    '''
    customized json serializer that handles numpy objects
    '''
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NpEncoder, self).default(obj)

if __name__ == '__main__':
    # set up the first parser to get the required epidemic parameter file
    p_params_file = argparse.ArgumentParser(add_help=False)
    p_params_file.add_argument("-pf", "--parameter-file", required=True, type=str, help="A JSON file with epidemic parameters. This can either be a abstract model file or complete model file. For abstract, you need also provide the epidemic rates")
    args, _ = p_params_file.parse_known_args()
    parameter_filename = args.parameter_file
    if not os.path.isfile(parameter_filename):
        print('Please provide a JSON file with epidemic parameter configurations as the first argument!')
        sys.exit(1)

    print('This file will be used: {}'.format(parameter_filename))
    epidemic_parameters = json.load(open(parameter_filename))

    tracing_parameters = json.load(open(DEFAULT_TRACING_PARAMS_FILE))
    population_parameters = json.load(open(DEFAULT_POPULATION_PARAMS_FILE))

    # now we use a second argparse to check other required parameters
    p = argparse.ArgumentParser()
    p.add_argument("-seed", type=int, default=-1, help='seed (if not specified seed is random')
    p.add_argument("-nt", "--network-type", type=str, default='er', help='specify the network simulator to use')
    p.add_argument("-np", "--network-params", type=str, default=None, help='specify the network parameter to use')
    p.add_argument("-oracle_type", type=str, default=None, help='type of oracle to use')
    p.add_argument("-num_tracers", type=int, default=0, help='number of tracers')
    p.add_argument("-nxIts", type=int, default=1, help='number of iterations of random networks')
    p.add_argument("-simIts", type=int, default=1, help='number of simulations per random network')
    p.add_argument("-of", type=str, default=None, help='path to the file to save results')
    p.add_argument("-return_full_data", type=str2bool, default=True, help='Determine if individual transimissions are tracked (This will introduce additional computational overhead)')
    p.add_argument("-parallel", type=str2bool, default=False, help='Run simulation in parallel.')
    p.add_argument("-num_nodes", type=int, default=0, help='number of first set of nodes selected for computing secondary infections')
    p.add_argument("-reintroduction_nums", type=int, default=0, help='number of infections at each re-introduction')
    for param in tracing_parameters:
        if isinstance(tracing_parameters[param], (int, float, bool)):
            p.add_argument('-'+param, type=type(tracing_parameters[param]), default=None)
    for param in population_parameters:
        if isinstance(population_parameters[param], (int, float, bool)):
            p.add_argument('-'+param, type=type(population_parameters[param]), default=None)

    # if abstract model file is provided, we need to get all rates
    if epidemic_parameters['type'] == 'abstract':
        print('An abstract model file is provided.')
        for param in epidemic_parameters['required_rates']:
            p.add_argument('-'+param, type=float, required=True, default=None)
    args, _ = p.parse_known_args()

    # regenerate parameters based on given parameters
    parameters = load_parameters(parameter_filename, DEFAULT_TRACING_PARAMS_FILE, DEFAULT_POPULATION_PARAMS_FILE, new_paramters=vars(args))
    print('The following parameters will be applied:\n{}'.format(parameters))

    if args.seed != -1:
        print('Using random seed: {}'.format(args.seed))
        random.seed(args.seed)
        np.random.seed(args.seed)
    # simulation function
    print('Start simulating loops...\n')
    infected_totals, max_infected, days2end, quarantined_totals,\
    average_secondary_infections, prob_positive_CT, prob_positive_RT,\
    ts, all_status_histories, full_transimissions, full_init_infections,\
    top_10_communities_infection, positive_per_day, tracers_per_day, tests_correlation_records,\
    network_sizes, list_secondary_infection_counts = simulation_loop(**parameters)

    print_results(parameters['N'], infected_totals, max_infected,
                  days2end, quarantined_totals, average_secondary_infections,
                  prob_positive_CT, prob_positive_RT, top_10_communities_infection, network_sizes)

    statuses = []
    for sta in ['I', 'E', 'Ho', 'Is', "Ia","Is", "Ps", "Ea","Es","Ho","ICU","H"]:
        if sta in all_status_histories[0]:
            statuses.append(sta)
    daily_counts = derive_daily_counts(all_status_histories, ts=ts, status_to_track=statuses)

    if args.of is not None:
        # if an output file path is provided, then we save the simulattion
        res = {
            "infected_totals": infected_totals,
            "max_infected": max_infected,
            "days2end": days2end,
            "quarantined_totals": quarantined_totals,
            "average_secondary_infections": average_secondary_infections,
            "prob_positive_CT": prob_positive_CT,
            "prob_positive_RT": prob_positive_RT,
            "positive_per_day": positive_per_day,
            "tracers_per_day": tracers_per_day,
            "daily_counts": daily_counts,
            "tests_correlation_records": tests_correlation_records,
            "top_10_communities_infection": top_10_communities_infection,
            "network_sizes": network_sizes
        }
        with open(args.of, 'w') as fp:
            json.dump(res, fp, cls=NpEncoder)
