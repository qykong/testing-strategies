import json
import random
import sys
from collections import defaultdict

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

from pyscripts.simulations import Gillespie_simple_contagion
from pyscripts.tracer import DefaultTracer

def generate_init_infections(nI, N, initial_infection_status):
    # Initial infection
    IC = defaultdict(lambda: 'S')

    for rn in random.sample(range(N), nI):
        IC[rn] = initial_infection_status

    return IC


def generate_contact_networks(N, network_type, network_params=None):
    if type(network_params) is str:
        network_params = json.loads(network_params)
    if network_type == 'er' or network_type == 'er_gamma':
        if network_params is not None:
            if 'p' in network_params:
                G = nx.fast_gnp_random_graph(N, network_params['p'])
            elif 'p_without_N' in network_params:
                G = nx.fast_gnp_random_graph(N, network_params['p_without_N']/(N-1))
        else:
            G = nx.fast_gnp_random_graph(N, 3. / (N - 1))
        G.remove_nodes_from(list(nx.isolates(G)))
    elif network_type == 'nbinomial':
        dist = pd.read_csv(network_params['dist_file'])
        degrees = np.random.choice(dist.degree.to_numpy(), p=dist.prob.to_numpy(), size=N)
        if sum(degrees) % 2 == 1:
            degrees[0] += 1
        G = nx.configuration_model(degrees)
        G = nx.Graph(G) # remove parallel edges
        G.remove_edges_from(nx.selfloop_edges(G)) # remove self loops
        print(nx.number_connected_components(G), ' connected component')
    elif network_type == 'aleta':
        G = nx.read_edgelist(open('data/aleta_network.csv', 'rb'))
        G = nx.convert_node_labels_to_integers(G)
    else:
        print('unknown network type!')
        sys.exit(1)
    return G


def load_parameters(file_path, tracing_params_file=None, population_params_file=None, new_paramters=None):
    '''
        read parameters from a json file and also set up the epi
        graphs for the EoN library
    '''
    if tracing_params_file is not None:
        with open(tracing_params_file) as f:
            tracing_params = json.load(f)
    if population_params_file is not None:
        with open(population_params_file) as f:
            population_params = json.load(f)
    with open(file_path) as f:
        parameters = json.load(f)
        if tracing_params_file is not None:
            parameters = {**tracing_params, **parameters}
        if population_params_file is not None:
            parameters = {**population_params, **parameters}
        if new_paramters is not None:
            for param in new_paramters:
                if new_paramters[param] is not None:
                    parameters[param] = new_paramters[param]
        # build epi graphs
        gH = nx.DiGraph()
        for item in parameters['gH']:
            if type(item) == list:
                gH.add_edge(item[0], item[1], rate=parameters[item[2]])
            elif type(item) == str:
                gH.add_node(item)
        J = nx.DiGraph()
        for item in parameters['J']:
            J.add_edge((item[0][0], item[0][1]), (item[1][0], item[1][1]), rate=parameters[item[2]])
        parameters['gH'] = gH
        parameters['J'] = J
        parameters['return_statuses'] = tuple(parameters['return_statuses'])
        
        return parameters

SAVED_WEIGHTS = {}
def process_gamma(G, J, network_params):
    # this functions is for gamma distributed infection rates
    if type(network_params) is str:
        network_params = json.loads(network_params)
    global SAVED_WEIGHTS
    if len(SAVED_WEIGHTS) > 0:
        SAVED_WEIGHTS = {}
    if not G.is_directed():
        G = G.to_directed()

    rate_IS2II = J.adj[('I','S')][('I','I')]['rate']
    def rf(G2, source, target, **nbr_kwargs):
        nonlocal rate_IS2II
        nonlocal network_params

        global SAVED_WEIGHTS
        if source not in SAVED_WEIGHTS:
            SAVED_WEIGHTS[source] = np.random.gamma(network_params['dispersion'], rate_IS2II/network_params['dispersion'], 1)[0]
        return max(SAVED_WEIGHTS[source], 10**(-7))
    
    J.adj[('I','S')][('I','I')]['rate'] = 1
    J.adj[('I','S')][('I','I')]['rate_function'] = rf
    return G, J

def simulation_loop(nxIts, simIts, nI, N, gH, J, return_statuses, maxdays, num_tracers, pC2T, pQ, p2bCTs,
                    p2bRTs, positive_instantaneous, positive_induced, infected, recovered, to_trace_instantaneous, pCT,
                    tlist_size, oracle=None, oracleTracer=None, return_full_data=True, network_type='er', network_params=None,
                    infected_to_count=['I'], initial_infection_status='I', oracle_type=None, parallel=False,
                    num_nodes=0, reintroduction_nums=0, p2bPITs=dict(), pRT=1, **kwargs):
    '''
    **return_full_data** determines if each infection event is recorded or not
    '''
    #### Vars for results ####
    infected_totals = []
    days2end = []
    max_infected = []
    quarantined_totals = []
    prob_positive_CT = []
    prob_positive_RT = []
    average_secondary_infections = []
    # I = []
    ts = []
    pds = []
    tds = []
    all_status_histories = []  # keep all status histories for plotting counts in chosen compartments
    full_transimissions = []
    full_init_infections = []
    top_10_communities_infection = []
    tests_correlation_records = []
    network_sizes = []
    list_secondary_infection_counts = []
    #### Simulations loop ####
    for its in range(nxIts):
        G = generate_contact_networks(N, network_type, network_params)
        # hack for gamma distribution
        if parallel:
            import multiprocessing
            from functools import partial
            if network_type == 'er_gamma':
                raise Exception('cannot run er_gamma in parallel!')
            func = partial(run_single_simulation, G, J, N, gH, infected, initial_infection_status, maxdays,
                                                                 nI, num_tracers, oracle, oracleTracer, oracle_type, p2bCTs,
                                                                 p2bRTs, pC2T, pCT, pQ, positive_induced,
                                                                 positive_instantaneous, recovered, return_full_data,
                                                                 return_statuses, tlist_size, to_trace_instantaneous,
                                                                 reintroduction_nums, p2bPITs, pRT)
            with multiprocessing.Pool(processes=40) as p:
                all_results = p.map(func, range(simIts))
        else:
            all_results = []
            for j in range(simIts):
                if network_type == 'er_gamma':
                    G, J = process_gamma(G, J, network_params)
                init_infections, results = run_single_simulation(G, J, N, gH, infected, initial_infection_status, maxdays,
                                                                 nI, num_tracers, oracle, oracleTracer, oracle_type, p2bCTs,
                                                                 p2bRTs, pC2T, pCT, pQ, positive_induced,
                                                                 positive_instantaneous, recovered, return_full_data,
                                                                 return_statuses, tlist_size, to_trace_instantaneous,
                                                                 reintroduction_nums, p2bPITs, pRT, j)
                all_results.append((init_infections, results))
                print("Done one simulation")
        for init_infections, results in all_results:
            # process results
            if return_full_data:
                t, td, pd, q, inf, tests_correlation_record, transmissions = results[:1] + results[-6:]
            else:
                t, td, pd, q, inf, tests_correlation_record = results[:1] + results[-5:]
            status_histories = {return_status: results[1 + i] for i, return_status in enumerate(return_statuses)}
            #### Save results ####
            infected_totals.append(len(inf))
            # get the time of last infection event
            index_of_last_infection_event = np.where(status_histories['S'] == np.min(status_histories['S']))[0][0]
            days2end.append(t[index_of_last_infection_event])
            infected_counts = sum(status_histories[status] for status in infected_to_count)
            max_infected.append(max(infected_counts))
            quarantined_totals.append(len(q))
            ts.append(t)
            pds.append(pd)
            tds.append(td)
            all_status_histories.append(status_histories)
            network_sizes.append(G.number_of_nodes())
            tests_correlation_records.append(tests_correlation_record)

            ppCT = 0.
            ppRT = 0.
            nppCT = 0
            nppRT = 0
            for k in range(len(td)):
                (cts, rts) = td[k]
                (pct, prt) = pd[k]
                if cts > 0:
                    ppCT += pct / float(cts)
                    nppCT += 1
                if rts > 0:
                    ppRT += prt / float(rts)
                    nppRT += 1
            if nppCT > 0:
                prob_positive_CT.append(ppCT / float(nppCT))
            if nppRT > 0:
                prob_positive_RT.append(ppRT / float(nppRT))
            # compute average secondary infections
            if return_full_data:
                # import ipdb; ipdb.set_trace()
                if num_nodes == -1:
                    # get secondary infections before peak
                    infections_c = np.diff(infected_counts)
                    infections_c[infections_c < 0] = 0
                    infections_c = np.concatenate([[0], infections_c])
                    daily_counts = np.unique(np.array(t // 1)[infections_c==1], return_counts = True)
                    
                    if len(daily_counts[1]) > 0:
                        max_day = daily_counts[0][np.argmax(daily_counts[1])]
                    else:
                        max_day = 0
                    secondary_infection_counts = {node: 0. for time, _, node in transmissions if time < max_day + 1}
                else:
                    secondary_infection_counts = {node: 0. for _, _, node in transmissions[:num_nodes]}

                for time, source, node in transmissions:
                    if source in secondary_infection_counts:
                        secondary_infection_counts[source] += 1
                list_secondary_infection_counts.append(list(secondary_infection_counts.values()))
                full_init_infections.append(init_infections)

            # check top 10 communities' infection status
            if not G.is_directed():
                Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
                GI_counts = []
                for k in range(min(10, len(Gcc))):
                    GI_counts.append(len(Gcc[k].intersection(inf))/len(Gcc[k]))
                top_10_communities_infection.append(GI_counts)

    return infected_totals, max_infected, days2end, quarantined_totals, \
    average_secondary_infections, prob_positive_CT, prob_positive_RT, ts, all_status_histories, full_transimissions, \
           full_init_infections, top_10_communities_infection, pds, tds, tests_correlation_records, network_sizes, list_secondary_infection_counts


def run_single_simulation(G, J, N, gH, infected, initial_infection_status, maxdays, nI, num_tracers, oracle,
                          oracleTracer, oracle_type, p2bCTs, p2bRTs, pC2T, pCT, pQ, positive_induced,
                          positive_instantaneous, recovered, return_full_data, return_statuses, tlist_size,
                          to_trace_instantaneous, reintroduction_nums, p2bPITs, pRT, j):
    IC = generate_init_infections(round(nI * G.number_of_nodes() / N), G.number_of_nodes(), initial_infection_status)
    init_infections = list(IC.keys())
    # keep oracleTracer and oracle for backward compatability
    # if oracle_type is given, ignore other variables
    if oracle_type is None:
        if oracle is None and oracleTracer is None:
            raise Exception('all oracle settings are None!')
        oracle_type = 'forward'
        if oracle:
            oracle_type = 'globaloracle'
        elif oracleTracer:
            oracle_type = 'oracletracer'
    if oracle_type == 'random':
        pCT = 0
    tracer = DefaultTracer(num_tracers, pC2T, pQ, p2bCTs, p2bRTs,
                           positive_instantaneous, positive_induced,
                           tlist_size, to_trace_instantaneous, oracle_type, pCT, p2bPITs, pRT)
    results = Gillespie_simple_contagion(G, gH, J, IC, return_statuses, 0, maxdays,
                                         infected, recovered, tracer=tracer, return_full_data=return_full_data, reintroduction_nums=reintroduction_nums)
    return init_infections, results


def print_results(N, infected_totals, max_infected, days2end, quarantined_totals, average_secondary_infections,
                  prob_positive_CT, prob_positive_RT, top_10_communities_infection, network_sizes):
    #### Print Results ####
    print("Average Ratio of total infected", np.mean(infected_totals / np.array(network_sizes)))
    print("Average Ratio of max infected", np.mean(max_infected / np.array(network_sizes)))
    print("Average Days to Epidemic End", np.mean(days2end))
    print("Average Ratio of quarantined", np.mean(quarantined_totals / np.array(network_sizes)))
    if len(average_secondary_infections) > 0:
        print("Average Secondary Infection Counts", np.mean(average_secondary_infections))
    else:
        print("No secondary infection counted")
    if len(prob_positive_CT) > 0:
        print("Average Ratio of Positive Contact Tracing", np.mean(prob_positive_CT))
    else:
        print("No Contact Tracing Performed")
    if len(prob_positive_RT) > 0:
        print("Average Ratio of Positive Random Tracing", np.mean(prob_positive_RT))
    else:
        print("No Random Tracing Performed")

def derive_daily_counts(all_status_histories, ts, status_to_track=['I']):
    daily_changes = {}
    for status in status_to_track:
        daily_changes[status] = []
        for t, status_histories in zip(ts, all_status_histories):
            floored_time = np.floor(t)
            max_day = np.max(floored_time)

            # compute daily changes of the chosen status
            status_history = status_histories[status]
            status_history_updates = np.concatenate([[0], np.diff(status_history)])
            daily = []
            for d in range(int(max_day)+1):
                daily.append((sum(status_history_updates[np.logical_and(floored_time == d,
                                                                    status_history_updates < 0)]),
                              sum(status_history_updates[np.logical_and(floored_time == d,
                                                                        status_history_updates > 0)])
                              ))
            daily_changes[status].append(daily)
    return daily_changes
