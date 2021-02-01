# the original code is from https://github.com/springer-math/Mathematics-of-Epidemics-on-Networks
# we modified the code to include tracing and testing strategies.
import math
import random
from collections import Counter, defaultdict

import EoN
import networkx as nx
import numpy as np

from pyscripts.listdict import _ListDict_


def Gillespie_simple_contagion(G, spontaneous_transition_graph,
                               nbr_induced_transition_graph, IC, return_statuses, tmin, tmax, infected, recovered,
                               tracer, spont_kwargs=None, nbr_kwargs=None, return_full_data=False, reintroduction_nums=0):
    r'''
    Performs simulations for epidemics, allowing more flexibility than SIR/SIS.

    This does not handle complex contagions.  It assumes that when an individual
    changes status either he/she has received a "transmission" from a *single*
    neighbor or he/she is changing status independently of any neighbors.  So
    this is like SIS or SIR.  Other examples would be like SEIR, SIRS, etc

    There is an example below demonstrating an SEIR epidemic.

    We allow for nodes to undergo two types of transitions.  They can be:

    - spontaneous - so a node of status A becomes B without any neighbor's
      influence

    - neighbor-induced - so an edge between status A and status B nodes suddenly
      becomes an edge between a status A and a status C node because the status
      B node changes status due to the existence of the edge.  (in principle,
      both could change status, but the code currently only allows the second
      node to change status).

    Both types of transitions can be represented by weighted directed graphs.
    We describe two weighted graphs whose nodes represent the possible statuses
    and whose edges represent possible transitions of statuses.

    - The spontaneous transitions can be represented by a graph whose nodes are
      the possible statuses and an edge from ``'A'`` to ``'B'`` represent that in the
      absence of any influence from others an indivdiual of status ``'A'``
      transitions to status ``'B'`` with default rate given by the weight of the edge
      in this spontaneous transition graph.  The rate may be modified by
      properties of the nodes in the contact network ``G``.

    - The neighbor-induced transitions can be represented by a "transitions
      graph" whose nodes are length-2 tuples.  The first entry represents the
      first individual of a partnership and the second represents the second
      individual.  **only the second individual changes status**.  An edge in
      the transitions graph from the node ``('A', 'B')`` to the node ``('A', 'C')``
      represents that an ``'AB'`` partnership in the contact network can cause the
      second individual to transition to status ``'C'``.  The weight of the edge
      in the represents the default transition rate.  The rate may be modified
      by properties of the nodes or the edge in the contact network ``G``.

    [for reference, if you look at Fig 4.3 on pg 122 of Kiss, Miller & Simon
    the graphs for **SIS** would be:
        ``spontaneous_transition_graph``: ``'I'``-> ``'S'`` with the edge weighted by ``gamma`` and

        ``nbr_induced_transition_graph``: ``('I', 'S')`` -> ``('I', 'I')`` with the edge weighted by ``tau``.

    For **SIR** they would be:
        ``spontaneous_transition_graph``: ``'I'``->``'R'`` with weight ``gamma`` and

        ``nbr_induced_transition_graph``: ``('I', 'S')`` -> ``('I', 'I')`` with rate ``tau``.
        ]

    These graphs must be defined and then input into the algorithm.

    It is possible to weight edges or nodes in the contact network ``G`` (that is,
    not the 2 directed networks defined above, but the original contact
    network) so that some of these transitions have different rates for
    different individuals/partnerships.  These are included as attributes in the
    contact network.  In the most general case, the transition rate depends
    on some function of the attributes of the nodes or edge in the contact
    network ``G``.

    There are two ways we introduce individual or pair-level heterogeneity in
    the population.  The first way is through introducing weights to individuals
    or their contacts which simply multiply the default rate.  The second is
    through including a function of nodes or pairs of nodes which will explicitly
    calculate the scaling factor, possibly including more information than
    we can include with the weights.  For a given transition, you can use at
    most one of these.

    - We first describe examples of weighting nodes/edges in the population

    So for the SIR case, if some people have higher recovery rate, we might
    define a node attribute ``'recovery_weight'`` for each node in ``G``, and the
    recovery would occur with rate G.node[node]['recovery_weight']*gamma.  So
    a value of 1.1 would correspond to a 10% increased recovery rate.  Since I
    don't know what name you might choose for the weight label as I write this
    algorithm, in defining the spontaneous transition graph (H), the ``'I'``->``'R'``
    edge would be given an attribute ``'weight_label'`` so that
    ``H.adj['I']['R']['weight_label'] = 'recovery_weight'``.
    If you define the attribute ``'weight_label'`` for an edge in ``H``, then it will be
    assumed that every node in ``G`` has a corresponding weight.  If no attribute
    is given, then it is assumed that all transitions happen with the original
    rate.

    We similarly define the weight_labels as edge attributes in the
    neighbor-induced transition graph.  The edges of the graph ``'G'`` have a
    corresponding ``G[u,v]['transmission_weight']``

    - Alternately we might introduce a function.

    So for the SIR case if the recovery rate depends on two attributes of a node
    (say, age and gender), we define a function ``rate_function(G,node)`` which
    will then look at G.node[node]['age'] and G.node[node]['gender'] and then
    return a factor which will be multiplied by ``gamma`` to give the recovery
    rate.  Similarly, if we are considering a neighbor induced transition and
    the rate depends on properties of both nodes we define another function
    ``rate_function(G,u,v)`` which may use attributes of u or v or the edge
    to find the appropriate scaling factor.



    :Arguments:

    **G** NetworkX Graph
        The underlying contact network.  If ``G`` is directed, we assume that
        "transmissions" can only go in the same direction as the edge.

    **spontaneous_transition_graph** Directed networkx graph
        The nodes of this graph are the possible statuses of a node in ``G``.
        An edge in this graph is a possible transition in ``G`` that occurs
        without any influence from neighbors.

        An edge in this directed graph is labelled with attributes

          - ``'rate'``   [a number, the default rate of the transition]
          - ``'weight_label'``  (optional) [a string, giving the label of
                            a **node** attribute in the contact network ``G``
                            that scales the transition rate]
          - ``'rate_function'`` (optional not combinable with ``'weight_label'``
                            for some edge.)
                            [a user-defined function of the contact network
                            and node that will scale the transition rate.
                            This cannot depend on the statuses of any nodes -
                            we must be able to calculate it once at the
                            beginning of the process.]  It will be called as
                            ``rate_function(G, node,**spont_kwargs)`` where
                            ``spont_kwargs`` is described below.

        Only one of ``'weight_label'`` and ``'rate_function'`` can be given.

        In the description below, let's use
          - ``rate = spontaneous_transition_graph.adj[Status1][Status2]['rate']``
          - ``weight_label = spontaneous_transition_graph.adj[Status1][Status2]['weight_label']``
          - ``rate_function = spontaneous_transition_graph.adj[Status1][Status2]['rate_function']``

        For a node ``u`` whose status is ``Status1``, the rate at which ``u``
        transitions to ``Status2`` is
          - ``rate``    if neither ``weight_label`` nor ``rate_function`` is defined.
          - ``rate*G.node[u][weight_label]`` if ``weight_label`` is defined.
          - ``rate*rate_function(G, u, **spont_kwargs)`` if ``rate_function`` is
             defined.


        So for example in the case of an SIR disease, this would be a graph
        with an isolated node ``'S'`` and an edge from node ``'I'`` to ``'R'`` with
        ``rate`` equal to ``gamma`` (the recovery rate).  It would not actually
        be necessary to have the node ``'S'``.


        Note that the rate_function version is more general, and in principle
        we don't need the weight_label option.  It's here for backwards
        compatibility and general simplicity purposes.


    **nbr_induced_transition_graph** Directed networkx graph

        The nodes of this graph are tuples with possible statuses of nodes
        at the end of an edge. The first node in the tuple is the node that
        could be affecting the second.  So for example for the SIR model
        we would expect a node ``('I', 'S')`` with an edge to ``('I', 'I')``.

        An edge in this directed graph is labelled with attributes

          - ``'rate'``   [a number, the default rate of the transition]
          - ``'weight_label'``  (optional) [a string, giving the label of
                            an *edge** attribute in the contact network ``G``
                            that scales the transition rate]
          - ``'rate_function'`` (optional not combinable with ``'weight_label'``
                            for some edge.)
                            [a user-defined function of the contact network
                            and source and target nodes that will scale the
                            transition rate.
                            This cannot depend on the statuses of any nodes -
                            we must be able to calculate it once at the
                            beginning of the process.
                            It will be called as
                            ``rate_function(G, source, target, *nbr_kwargs)``]

        Only one of ``'weight_label'`` and ``'rate_function'`` can be given.

        In the description below, let's use

         - ``rate = spontaneous_transition_graph.adj[Status1][Status2]['rate']``
         - ``weight_label = spontaneous_transition_graph.adj[Status1][Status2]['weight_label']``
         - ``rate_function = spontaneous_transition_graph.adj[Status1][Status2]['rate_function']``

        If the transition is (A,B) to (A,C) and we have a source node of
        status A joined to a target node of status B then the target node
        transitions from B to C with rate
          - ``rate``    if neither ``weight_label`` nor ``rate_function`` is defined.
          - ``rate*G.adj[source][target][weight_label]`` if ``weight_label`` is defined.
          - ``rate*rate_function(G, source, target, **nbr_kwargs)`` if ``rate_function`` is defined

        So for example in the case of an SIR disease with transmission rate
        tau, this would be a graph with an edge from the node ``('I','S')`` to
        ``('I', 'I')``.  The attribute ``'rate'`` for the edge would be ``tau``.

    **IC** dict
        states the initial status of each node in the network.

    **return_statuses** list or other iterable (but not a generator)
        The statuses that we will return information for, in the order
        we will return them.

    **tmin** number (default 0)
        starting time

    **tmax** number (default 100)
        stop time

    **spont_kwargs**  dict or None (default None)
        Any parameters which might be needed if the user has defined a rate
        function for a spontaneous transition.
        If any of the spontaneous transition rate functions accepts these,
        they all need to (even if not used).  It's easiest to define the
        function as def f(..., **kwargs)

    **nbr_kwargs** dict or None (default None)
        Any parameters which might be needed if the user has defined a rate
        function for a neighbor-induced transition.
        If any of the neighbor-induced transition rate functions accepts these,
        they all need to (even if not used).  It's easiest to define the
        function as def f(..., **kwargs)

    **return_full_data** boolean
        Tells whether to return a Simulation_Investigation object or not

    **sim_kwargs** keyword arguments
        Any keyword arguments to be sent to the Simulation_Investigation object
        Only relevant if ``return_full_data=True``

    :Returns:

    **(times, status1, status2, ...)**  tuple of numpy arrays
        first entry is the times at which events happen.
        second (etc) entry is an array with the same number of entries as ``times``
        giving the number of nodes of status ordered as they are in ``return_statuses``


    :SAMPLE USE:

    This does an SEIR epidemic.  It treats the nodes and edges as weighted.
    So some individuals have a higher E->I transition rate than others, and some
    edges have a higher transmission rate than others.  The recovery rate
    is taken to be the same for all nodes.

    There are more examples in the
    online documentation at :ref:``simple-contagion-section``.

    ::

        import EoN
        import networkx as nx
        from collections import defaultdict
        import matplotlib.pyplot as plt
        import random

        N = 100000
        G = nx.fast_gnp_random_graph(N, 5./(N-1))

        #they will vary in the rate of leaving exposed class.
        #and edges will vary in transition rate.
        #there is no variation in recovery rate.

        node_attribute_dict = {node: 0.5+random.random() for node in G.nodes()}
        edge_attribute_dict = {edge: 0.5+random.random() for edge in G.edges()}

        nx.set_node_attributes(G, values=node_attribute_dict, name='expose2infect_weight')
        nx.set_edge_attributes(G, values=edge_attribute_dict, name='transmission_weight')


        H = nx.DiGraph()
        H.add_node('S')
        H.add_edge('E', 'I', rate = 0.6, weight_label='expose2infect_weight')
        H.add_edge('I', 'R', rate = 0.1)

        J = nx.DiGraph()
        J.add_edge(('I', 'S'), ('I', 'E'), rate = 0.1, weight_label='transmission_weight')
        IC = defaultdict(lambda: 'S')
        for node in range(200):
            IC[node] = 'I'

        return_statuses = ('S', 'E', 'I', 'R')

        t, S, E, I, R = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses,
                                                tmax = float('Inf'))

        plt.semilogy(t, S, label = 'Susceptible')
        plt.semilogy(t, E, label = 'Exposed')
        plt.semilogy(t, I, label = 'Infected')
        plt.semilogy(t, R, label = 'Recovered')
        plt.legend()

        plt.savefig('SEIR.png')
'''
    if spont_kwargs is None:
        spont_kwargs = {}

    if nbr_kwargs is None:
        nbr_kwargs = {}

    # re-introduction
    no_waves = 0
    # infected total individuals
    infected_total = set()
    #######################

    # dict node->status
    status = {node: IC[node] for node in G.nodes()}

    if return_full_data:
        transmissions = []
        # dict node->(list of times, list of status)
        node_history = {node: ([tmin], [status[node]]) for node in G.nodes()}

    times = [tmin]
    data = {}
    # Counter of number of nodes per status, dict = status - node
    C = Counter(status.values())
    for return_status in return_statuses:
        data[return_status] = [C[return_status]]  # list of Counter - initialized here
        
    # for measure daily tests correlation
    tests_correlation_record = {'in_positive_induced_daily': [],
                                'tested_positive_total_daily_start_of_day': [],
                                'quarantined_total_daily_start_of_day': [], # different from tested_positive_total when pQ=1
                                'tested_positive_daily_end_of_day': [],
                                'in_positive_instantaneous_daily': [],
                                'in_quarantine_daily_start_of_day': [],
                                'in_ctd_daily_start_of_day': []}

    # List of edges of spontaneous transitions
    spontaneous_transitions = list(spontaneous_transition_graph.edges())
    # List of edges of induced transitions ((A,B),(A,C))
    induced_transitions = list(nbr_induced_transition_graph.edges())
    potential_transitions = {}
    rate = {}  # intrinsic rate of a transition. Dict transition -> rate
    # weight_sum = defaultdict(lambda: 0)
    # weights = defaultdict(lambda: None)
    # max_weight = defaultdict(lambda: 0)
    get_weight = defaultdict(lambda: defaultdict(
        lambda: -1))  # weight increment defaults to -1 as ListDict is compiled with cython and needs to have a type

    # SET UP THE POSSIBLE EVENTS, STARTING WITH SPONTANEOUS
    for transition in spontaneous_transitions:
        rate[transition] = spontaneous_transition_graph.adj[transition[0]][transition[1]]['rate']
        if 'weight_label' in spontaneous_transition_graph.adj[transition[0]][transition[1]]:
            if 'rate_function' in spontaneous_transition_graph.adj[transition[0]][transition[1]]:
                raise EoN.EoNError('cannot define both "weight_label" and "rate_function" in',
                                   'spontaneous_transitions graph')
            wl = spontaneous_transition_graph.adj[transition[0]][transition[1]]['weight_label']
            get_weight[transition] = nx.get_node_attributes(G, wl)  # This is a dict mapping node to its weight.
            potential_transitions[transition] = _ListDict_(
                weighted=True)  # max_weight[transition] = max(get_weight[transition].values())
        elif 'rate_function' in spontaneous_transition_graph.adj[transition[0]][transition[1]]:
            rf = spontaneous_transition_graph.adj[transition[0]][transition[1]]['rate_function']
            get_weight[transition] = {node: rf(G, node, **spont_kwargs) for node in
                                      G}  # This is a dict mapping node to its weight.
            potential_transitions[transition] = _ListDict_(weighted=True)
        else:
            potential_transitions[transition] = _ListDict_()  # max_weight[transition]=1

    # print(spontaneous_transitions)
    # print([rate[transition] for transition in spontaneous_transitions])
    # CONTINUING SETTING UP POSSIBLE EVENTS, NOW WITH INDUCED TRANSITIONS.
    for transition in induced_transitions:
        if transition[0][0] != transition[1][0]:
            raise EoN.EoNError("transition {} -> {} not allowed: first node must keep same status".format(transition[0],
                                                                                                          transition[
                                                                                                              1]))
        rate[transition] = nbr_induced_transition_graph.adj[transition[0]][transition[1]]['rate']

        if 'weight_label' in nbr_induced_transition_graph.adj[transition[0]][transition[1]]:
            if 'rate_function' in nbr_induced_transition_graph.adj[transition[0]][transition[1]]:
                raise EoN.EoNError('cannot define both "weight_label" and "rate_function" in',
                                   'nbr_induced_transitions graph')
            wl = nbr_induced_transition_graph.adj[transition[0]][transition[1]]['weight_label']
            get_weight[transition] = nx.get_edge_attributes(G, wl)  # a dict mapping edge to its weight.
            # note if G is directed, then this has all edges.  If undirected
            # edge will appear only once.  But we may care about the opposite direction.
            # so we need to add those now.
            if not nx.is_directed(G):
                get_weight[transition].update(
                    {(source, target): G.adj[source][target][wl] for target, source in get_weight[transition]})
            potential_transitions[transition] = _ListDict_(weighted=True)
        elif 'rate_function' in nbr_induced_transition_graph.adj[transition[0]][transition[1]]:
            rf = nbr_induced_transition_graph.adj[transition[0]][transition[1]]['rate_function']
            edges = list(G.edges())

            get_weight[transition] = {(source, target): rf(G, source, target, **nbr_kwargs) for source, target in edges}
            if not nx.is_directed(G):
                get_weight[transition].update(
                    {(source, target): rf(G, source, target, **nbr_kwargs) for target, source in edges})
            potential_transitions[transition] = _ListDict_(weighted=True)
        else:
            potential_transitions[transition] = _ListDict_()

    # initialize all potential events to start.
    for node in G.nodes():
        if spontaneous_transition_graph.has_node(
                status[node]):  # and spontaneous_transition_graph.degree(status[node])>0:
            for transition in spontaneous_transition_graph.edges(status[node]):
                potential_transitions[transition].update(node, weight_increment=get_weight[transition][node])
                # weight increment defaults to None if not present

        for nbr in G.neighbors(node):
            # print(status[node],status[nbr])
            if nbr_induced_transition_graph.has_node((status[node], status[
                nbr])):  # and nbr_induced_transition_graph.degree((status[node],status[nbr])) >0:
                for transition in nbr_induced_transition_graph.edges((status[node], status[nbr])):
                    potential_transitions[transition].update((node, nbr),
                                                             weight_increment=get_weight[transition][(node, nbr)])
        if status[node] in infected:
            infected_total.add(node)
    t = tmin
    
    tracer.register_simulation_variables(G, spontaneous_transitions, induced_transitions,
                                         get_weight, status, potential_transitions, tmin, tmax)
    tracer.before_simulation()

    # NOW WE'RE READY TO GET STARTED WITH THE SIMULATING

    total_rate = sum(rate[transition] * potential_transitions[transition].total_weight() for transition in
                     spontaneous_transitions + induced_transitions)
    if total_rate > 0:
        delay = random.expovariate(total_rate)
    else:
        delay = float('Inf')
    t = t + delay
    ti = 0
    while total_rate > 0 and t < tmax:
        ### If new day, do tracing first ###
        if math.floor(t) >= tracer.tday + 1:
            tests_correlation_record['in_quarantine_daily_start_of_day'].append(len(tracer.quarantined))
            tests_correlation_record['in_ctd_daily_start_of_day'].append(len(tracer.ctd))
            tests_correlation_record['quarantined_total_daily_start_of_day'].append(len(tracer.quarantined_total))
            tests_correlation_record['tested_positive_total_daily_start_of_day'].append(len(tracer.tested_positive_total))
            tracer.action_on_new_day()
            # record daily in infected and tested positive cases
            tests_correlation_record['tested_positive_daily_end_of_day'].append(tracer.positives_day[-1])
            tests_correlation_record['in_positive_induced_daily'].append(sum(data[state][-1] for state in tracer.positive_induced))
            tests_correlation_record['in_positive_instantaneous_daily'].append(sum(data[state][-1] for state in tracer.positive_instantaneous))
            
            total_rate = sum(rate[transition] * potential_transitions[transition].total_weight() for transition in
                             spontaneous_transitions + induced_transitions)
            ###################################
        ####################################
        if math.floor(t) >= tracer.tday + 1:
            # this means we need to take into account another day of tracing resource
            continue
        times.append(t)
        ti += 1
        r = random.random()
        # Crappy thing should disappear
        # saved_transition = transition
        for transition in spontaneous_transitions + induced_transitions:
            if potential_transitions[transition].total_weight() > 0:
                saved_transition = transition
            r -= rate[transition] * potential_transitions[transition].total_weight() / total_rate
            if r < 0:
                break
        # if r > 0: # this should not happen
        #    transition = saved_transition
        # either node doing spontaneous or edge doing an induced event
        spontaneous = False
        if transition in spontaneous_transitions:
            spontaneous = True

        actor = potential_transitions[transition].choose_random()
        if spontaneous:
            modified_node = actor
            old_status = transition[0]
            new_status = transition[1]
            # node changes status
            n2t = modified_node
            if new_status in infected:
                infected_total.add(n2t)
            if new_status in recovered:
                tracer.node_recover(n2t)

            tracer.action_on_new_event(new_status, n2t)
            #############

        else:

            source, target = actor

            if source in tracer.quarantined or target in tracer.quarantined:
                print("This should not happen")

            modified_node = target
            n2t = modified_node
            old_status = transition[0][1]
            new_status = transition[1][1]
            if new_status in infected:
                infected_total.add(n2t)

            tracer.action_on_new_event(new_status, n2t)
            #############

            if return_full_data:
                transmissions.append((t, source, modified_node))

            # modified_node changes status

        status[modified_node] = new_status
        if return_full_data:
            node_history[modified_node][0].append(t)
            node_history[modified_node][1].append(new_status)

        # it might look like there is a cleaner way to do this, but what if it
        # happens that old_status == status[modified_node]???  This way still works.
        for x in data.keys():
            data[x].append(data[x][-1])
        if old_status in return_statuses:
            data[old_status][-1] -= 1
        if status[modified_node] in return_statuses:
            data[status[modified_node]][-1] += 1

        # UPDATE OUR POTENTIAL TRANSITIONS

        for transition in spontaneous_transitions:  # can probably make more efficient, but there aren't many
            # remove modified_node from any spontaneous lists
            # add modified_node to any spontaneous lists
            if transition[0] == old_status:
                potential_transitions[transition].remove(modified_node)
            if transition[0] == status[modified_node]:
                potential_transitions[transition].update(modified_node,
                                                         weight_increment=get_weight[transition][modified_node])
            # roundoff error can kill the calculation, but it's slow to do this right.
            # so we'll only deal with it if the value is small enough that roundoff
            # error might matter.
            if potential_transitions[transition].weighted and potential_transitions[transition].total_weight() < 10 ** (
                    -7) and potential_transitions[
                transition].total_weight() != 0:
                potential_transitions[transition].update_total_weight()

        for transition in induced_transitions:
            if G.is_directed():
                for nbr in G.neighbors(modified_node):
                    # remove edge from any induced lists
                    # add edge to any induced lists

                    nbr_status = status[nbr]

                    if (modified_node, nbr) not in get_weight[transition]:
                        get_weight[transition][(modified_node, nbr)] = get_weight[transition][(nbr, modified_node)]
                    if transition[0] == (old_status, nbr_status):
                        potential_transitions[transition].remove((modified_node, nbr))
                    if modified_node not in tracer.quarantined and nbr not in tracer.quarantined:
                        if transition[0] == (status[modified_node], nbr_status):
                            potential_transitions[transition].update((modified_node, nbr),
                                                                     weight_increment=get_weight[transition][
                                                                         (modified_node, nbr)])
                for pred in G.predecessors(modified_node):
                    # remove edge from any induced lists
                    # add edge to any induced lists

                    pred_status = status[pred]

                    if (pred, modified_node) not in get_weight[transition]:
                        get_weight[transition][(pred, modified_node)] = get_weight[transition][(pred, modified_node)]
                    if transition[0] == (pred_status, old_status):
                        potential_transitions[transition].remove((pred, modified_node))
                    if modified_node not in tracer.quarantined and pred not in tracer.quarantined:
                        if transition[0] == (pred_status, status[modified_node]):
                            potential_transitions[transition].update((pred, modified_node),
                                                                     weight_increment=get_weight[transition][
                                                                         (pred, modified_node)])
            else:
                for nbr in G.neighbors(modified_node):
                    # remove edge from any induced lists
                    # add edge to any induced lists
                    nbr_status = status[nbr]

                    if (modified_node, nbr) not in get_weight[transition]:
                        get_weight[transition][(modified_node, nbr)] = get_weight[transition][(nbr, modified_node)]
                    elif (nbr, modified_node) not in get_weight[transition]:
                        get_weight[transition][(nbr, modified_node)] = get_weight[transition][(modified_node, nbr)]

                    if transition[0] == (nbr_status, old_status):
                        potential_transitions[transition].remove((nbr, modified_node))
                    if transition[0] == (old_status, nbr_status):
                        potential_transitions[transition].remove((modified_node, nbr))

                    if modified_node not in tracer.quarantined and nbr not in tracer.quarantined:
                        if transition[0] == (nbr_status, status[modified_node]):
                            potential_transitions[transition].update((nbr, modified_node),
                                                                     weight_increment=get_weight[transition][
                                                                         (nbr, modified_node)])
                        if transition[0] == (status[modified_node], nbr_status):
                            potential_transitions[transition].update((modified_node, nbr),
                                                                     weight_increment=get_weight[transition][
                                                                         (modified_node, nbr)])

            # roundoff error can kill the calculation, but it's slow to do this right.
            # so we'll only deal with it if the value is small enough that roundoff
            # error might matter.
            if potential_transitions[transition].weighted and potential_transitions[transition].total_weight() < 10 ** (-7) and potential_transitions[
                transition].total_weight() != 0: # 
                potential_transitions[transition].update_total_weight()

                
        # use for loop instead of comprehension here to speed up
        total_rate = 0
        for transition in spontaneous_transitions + induced_transitions:
            total_rate += rate[transition] * potential_transitions[transition].total_weight()

        # check here if re-introduction is enabled
        # hard coded number of waves
        if reintroduction_nums > 0 and total_rate == 0 and no_waves <= 4:
            transition = ('E', 'I')
            if transition not in potential_transitions:
                raise Exception('not implemented for SIR yet')
            # let's hard code the status for COVID and SARS
            s_nodes = []
            for node in G.nodes():
                if status[node] == 'S':
                    s_nodes.append(node)
            
            if len(s_nodes) == 0:
                reintroduction_nums = 0
            else:
                print('another wave:', no_waves, ' infected_total:', len(infected_total))
                selected_nodes = random.sample(s_nodes, min(reintroduction_nums, len(s_nodes)))
                for selected_node in selected_nodes:
                    status[selected_node] = 'E'
                    potential_transitions[transition].update(selected_node,
                                                         weight_increment=get_weight[transition][selected_node])
                no_waves += 1
                # roundoff error can kill the calculation, but it's slow to do this right.
                # so we'll only deal with it if the value is small enough that roundoff
                # error might matter.
                if potential_transitions[transition].weighted and potential_transitions[transition].total_weight() < 10 ** (
                        -7) and potential_transitions[
                    transition].total_weight() != 0:
                    potential_transitions[transition].update_total_weight()
                total_rate = rate[transition] * potential_transitions[transition].total_weight()

        if total_rate > 0:
            delay = random.expovariate(total_rate)
        else:
            delay = float('Inf')
        t += delay
    
    # plug no_wave here... as it's the easiest to do
    tests_correlation_record['no_waves'] = no_waves
    returnval = []
    times = np.array(times)
    returnval.append(times)
    for return_status in return_statuses:
        data[return_status] = np.array(data[return_status])
        returnval.append(data[return_status])
    returnval.append(tracer.tracers_day)
    returnval.append(tracer.positives_day)
    returnval.append(tracer.quarantined_total)
    returnval.append(infected_total)
    returnval.append(tests_correlation_record)
    if return_full_data:
        returnval.append(transmissions)
    return returnval
