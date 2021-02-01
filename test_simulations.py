# this script hosts some test cases to make sure all changes don't break simulations
import datetime

from pyscripts.utils import *

DEFAULT_TRACING_PARAMS_FILE = 'parameters/default_tracing_parameters.json'
DEFAULT_POPULATION_PARAMS_FILE = 'parameters/default_population_parameters.json'

random.seed(120)
N = 10000
tracing_parameters = {
    "oracle_type": 'none',
    "num_tracers": 0,
    "pC2T": 1,
    "pQ": 1,
    "pCT": 1,
    "tlist_size": 0,
    "nI": 20,
    "maxdays": 60,
    "nxIts": 1,
    "simIts": 30,
    "N": N,
    "infected_to_count": ["Is", "Ia"],
    "network_type": "poisson",
    "network_params": {"p": 4/(N-1)}
}

def single_test():
    begin = datetime.datetime.now()
    print(f'test COVID model: oracle_type={tracing_parameters["oracle_type"]}',
          f'N={tracing_parameters["N"]} num_tracers={tracing_parameters["num_tracers"]}',
          f'tlist={tracing_parameters["tlist_size"]}')
    parameters = load_parameters('parameters/complete/Aleta2020_covid19_model.json', new_paramters=tracing_parameters)
    infected_totals, max_infected, days2end, quarantined_totals,\
    average_secondary_infections, prob_positive_CT, prob_positive_RT,\
    ts, all_status_histories, full_transimissions, full_init_infections,\
    top_10_communities_infection, positive_per_day, tracers_per_day, tests_correlation_records,\
    network_sizes, list_secondary_infection_counts = simulation_loop(**parameters)
    print_results(parameters['N'], infected_totals, max_infected, days2end, quarantined_totals,
                  average_secondary_infections, prob_positive_CT, prob_positive_RT, top_10_communities_infection)
    end = datetime.datetime.now()
    print('This takes {} seconds'.format((end - begin).seconds))
#%%
single_test()

#%%
tracing_parameters['num_tracers'] = 100
single_test()

#%%
tracing_parameters['oracle_type'] = 'noneback'
single_test()

#%%
tracing_parameters['num_tracers'] = 1000
single_test()

#%%
tracing_parameters["oracle_type"] = 'oracleTracer'
tracing_parameters['num_tracers'] = 100
single_test()

#%%
tracing_parameters['tlist_size'] = 1
single_test()

#%%
tracing_parameters['num_tracers'] = 1000
single_test()

#%%
tracing_parameters["oracle_type"] = 'oracle'
tracing_parameters['num_tracers'] = 1000
tracing_parameters['tlist_size'] = 1
single_test()
