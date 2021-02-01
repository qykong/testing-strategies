import math
import random
import networkx as nx
from collections import OrderedDict
import numpy as np

class Tracer(object):
    """
    This abstract class defines required methods for a tracing method to be applied in the
    simulation
    """

    def __init__(self, num_tracers, pC2T, pQ, p2bCTs, p2bRTs, positive_instantaneous, positive_induced,
                 tlist_size, to_trace_instantaneous, p2bPITs=dict()):
        """
        **num_tracers** number of tracers available
        **positive_instantaneous** list of statuses that are instantaneous positives
        **positive_induced** list of statuses that become positive when traced
        **pCT** probability of tracer doing contact tracing
        **pC2T** probability of real contacts being discovered
        **pQ** probability of positive following quarantine recomendations
        **p2bCTs** dict of probabilities per status of being randomly traced (status->prob)
        **p2bRTs** dict of probabilities per status of being randomly traced (status->prob)
        **tlist_size** size of the tabu list that prevents nodes from being retraced
        """
        self.num_tracers = num_tracers
        self.pC2T = pC2T
        self.pQ = pQ
        self.p2bCTs = p2bCTs
        self.p2bRTs = p2bRTs
        self.p2bPITs = p2bPITs
        self.positive_instantaneous = positive_instantaneous
        self.positive_induced = positive_induced
        self.tlist_size = tlist_size
        self.to_trace_instantaneous = to_trace_instantaneous

        ###Tracing variables###
        # list of nodes to trace
        self.c2t = OrderedDict()

        # list of nodes traced
        self.ctd = set()
        # nodes that are quarantined
        self.quarantined = set()
        self.quarantined_total = set()
        self.tested_positive_total = set() # this is the same as quarantined_total if pQ=1

        # list of positives discovered per day (Ctpos,RTpos)
        self.positives_day = []

    def register_simulation_variables(self, G, spontaneous_transitions, induced_transitions, get_weight,
                                      status, potential_transitions, tmin, tmax):
        self.G = G
        self.spontaneous_transitions = spontaneous_transitions
        self.induced_transitions = induced_transitions
        self.get_weight = get_weight
        self.potential_transitions = potential_transitions
        self.status = status
        self.tmin = tmin
        self.tmax = tmax
        # tracing_day
        self.tday = tmin
        
        # save k-core values beforehand
        self.G_kcores = nx.core_number(G)

        # tabu list for tracing
        self.tlist = set()  # a set of nodes who cannot be traced due to tabu
        self.tlist_remove_time = {day: [] for day in
                                  range(
                                      tmax + 1)}  # key is the day to be removed and value is a list of nodes to remove

    def quarantine(self, n2t):
        # try to quarantine
        rt = random.random()
        if rt < self.pQ:  # the node follows quarantine suggestion
            self.quarantined.add(n2t)
            if n2t not in self.quarantined_total:
                self.quarantined_total.add(n2t)
            # update potential transitions
            old_status = self.status[n2t]
            for transition in self.induced_transitions:
                if self.G.is_directed():
                    for nbr in self.G.neighbors(n2t):
                        # remove edge from any induced lists
                        nbr_status = self.status[nbr]
                        if (n2t, nbr) not in self.get_weight[transition]:
                            self.get_weight[transition][(n2t, nbr)] = self.get_weight[transition][(nbr, n2t)]
                        if transition[0] == (old_status, nbr_status):
                            self.potential_transitions[transition].remove((n2t, nbr))
                    for pred in self.G.predecessors(n2t):
                        # remove edge from any induced lists
                        pred_status = self.status[pred]
                        if (pred, n2t) not in self.get_weight[transition]:
                            self.get_weight[transition][(pred, n2t)] = self.get_weight[transition][(pred, n2t)]
                        if transition[0] == (pred_status, old_status):
                            self.potential_transitions[transition].remove((pred, n2t))
                else:
                    for nbr in self.G.neighbors(n2t):
                        # remove edge from any induced lists
                        nbr_status = self.status[nbr]
                        if (n2t, nbr) not in self.get_weight[transition]:
                            self.get_weight[transition][(n2t, nbr)] = self.get_weight[transition][(nbr, n2t)]
                        elif (nbr, n2t) not in self.get_weight[transition]:
                            self.get_weight[transition][(nbr, n2t)] = self.get_weight[transition][(n2t, nbr)]

                        if transition[0] == (nbr_status, old_status):
                            self.potential_transitions[transition].remove((nbr, n2t))
                        if transition[0] == (old_status, nbr_status):
                            self.potential_transitions[transition].remove((n2t, nbr))

                # roundoff error can kill the calculation, but it's slow to do this right.
                # so we'll only deal with it if the value is small enough that roundoff
                # error might matter.
                if self.potential_transitions[transition].weighted and self.potential_transitions[
                    transition].total_weight() < 10 ** (
                        -7) and self.potential_transitions[
                    transition].total_weight() != 0:
                    self.potential_transitions[transition].update_total_weight()

    def node_recover(self, node):
        if node in self.ctd:
            self.ctd.remove(node)
        if node in self.quarantined:
            self.quarantined.remove(node)

    def before_simulation(self):
        """
        hook on actions before the simulation starts
        :return:
        """
        pass

    def action_on_new_day(self):
        """
        do whatever is required on the new simulation day to trace over the disease spreading network.
        :return:
        """
        pass

    def action_on_new_event(self):
        """
        this function is mainly for oracle types where you can do something whenever there's a new event happens
        :return: None
        """
        pass


class DefaultTracer(Tracer):

    def __init__(self, num_tracers, pC2T, pQ, p2bCTs, p2bRTs, positive_instantaneous, positive_induced,
                 tlist_size, to_trace_instantaneous, oracle_type, pCT, p2bPITs=dict(), pRT=1, trace_depth=1):
        """
        :param oracle_type: type of oracles to use, one of the ('none', 'oracle', 'oracleTracer')
        :param pCT: percentages of tracers doing contact tracing
        """
        super().__init__(num_tracers, pC2T, pQ, p2bCTs, p2bRTs, positive_instantaneous, positive_induced,
                         tlist_size, to_trace_instantaneous, p2bPITs)
        if oracle_type.lower() not in ['noneback', 'none', 'oracletracer', 'oracle', 'random', 'backoracletracer', 'kcorecontact']:
            raise Exception('oracle type not found!')

        # add a new variable for checking backward tracing depth
        self.node_trace_depth = {}
        self.oracle_type = oracle_type
        self.pCT = pCT
        self.pRT = pRT
        self.trace_depth = trace_depth
        
        # list of tracers per day (nCT,nRT)
        self.tracers_day = []

    def c2t_if_oracle(self, n2t):
        if self.oracle_type.lower() == 'oracle' and self.status in self.to_trace_instantaneous:
            if n2t not in self.ctd and n2t not in self.c2t and n2t not in self.tlist:
                rt = random.random()
                if rt < self.p2bCTs[self.status]:
                    self.c2t[n2t] = 1

    def find_out_contacts(self, n2t):
        if self.oracle_type.lower() != 'oracle':
            for nbr in self.G.neighbors(n2t):
                rt = random.random()
                if rt < self.pC2T:  # the neighbor node is discovered
                    if nbr not in self.ctd and nbr not in self.c2t and nbr not in self.tlist:
                        self.node_trace_depth[nbr] = self.node_trace_depth[n2t] + 1
                        if (self.oracle_type.lower() == 'noneback' or self.oracle_type.lower() == 'backoracletracer') and self.node_trace_depth[n2t] < self.trace_depth:
                            # backward tracing
                            self.c2t[nbr] = 1
                            self.c2t.move_to_end(nbr, last=False)
                        else:
                            # forward tracing
                            self.c2t[nbr] = 1

    def status_in_positive_instantaneous(self, new_status, n2t):
        if new_status in self.positive_instantaneous and \
                ((new_status in self.p2bPITs and random.random() < self.p2bPITs[new_status]) or \
                (new_status not in self.p2bPITs)):
            # n2t might be in ctd already as, for example, it moves from a compartment that tests possitive (Is) to
            # an positive instantaneous status, Ho
            if n2t not in self.ctd:
                self.ctd.add(n2t)
                if n2t not in self.node_trace_depth:
                    self.node_trace_depth[n2t] = 0
                self.find_out_contacts(n2t)
                self.quarantine(n2t)

    def orderC2T(self):
        pc2t = 0
        nodes = list(self.c2t.keys())
        for node in reversed(nodes):
            if self.status[node] in self.positive_induced:
                self.c2t.move_to_end(node, last=False)
                pc2t += 1
        return pc2t

    def before_simulation(self):
        ####Some nodes might already be positive
        for n2t in self.G.nodes():
            # if oracle then put "to_trace_instantaneous" people in c2t list
            self.c2t_if_oracle(n2t)
            self.status_in_positive_instantaneous(self.status[n2t], n2t)

    def action_on_new_day(self):
        self.tday = self.tday + 1
        if len(self.tlist_remove_time[self.tday]) > 0:
            self.tlist = self.tlist - set(self.tlist_remove_time[self.tday])
        # print(tday,len(c2t))
        ###### Start Contact/Random tracing ######
        CTpos = 0
        RTpos = 0
        if self.oracle_type.lower() == 'oracletracer' or self.oracle_type.lower() == 'backoracletracer':  # oracle tracing
            pc2t = self.orderC2T()
        else:
            pc2t = len(self.c2t)
        nCT = math.floor(self.pCT * self.num_tracers)
        nCT = min(pc2t, nCT)  # Best we can do, we only do CT to positives
        nRT = math.floor((self.num_tracers - nCT) * self.pRT)
        # We do a for here because we assume contacts take up tracing resource even if they do not accept testing
        for i in range(nCT):
            n2t = self.c2t.popitem(last=False)[0]
            if self.tlist_size > 0 and self.tday + self.tlist_size <= self.tmax:
                self.tlist.add(n2t)
                self.tlist_remove_time[self.tday + self.tlist_size].append(n2t)
            # try to trace it
            rt = random.random()
            if rt < self.p2bCTs[self.status[n2t]]:
                if self.status[n2t] in self.positive_induced:
                    CTpos += 1
                    self.tested_positive_total.add(n2t)
                    if n2t not in self.ctd:
                        self.ctd.add(n2t)
                        if n2t not in self.node_trace_depth:
                            self.node_trace_depth[n2t] = 0
                    else:
                        print("Also not good")
                    # find out contacts
                    self.find_out_contacts(n2t)
                    self.quarantine(n2t)
                    #############

        # It could be that there is not enough people to trace, we assume here a recovery/death model
        available_nodes_to_trace = set(self.G.nodes()) - self.ctd - self.tlist
        nRT = min(nRT, len(available_nodes_to_trace))
        # update tracers per day
        self.tracers_day.append((nCT, nRT))
        # do random tracing here
        if self.oracle_type.lower() == 'kcorecontact':
            available_nodes_to_trace_list = list(available_nodes_to_trace)
            available_nodes_to_trace_list_p = np.array([self.G_kcores[node] for node in available_nodes_to_trace_list])
            sampled_nodes = np.random.choice(available_nodes_to_trace_list, size=nRT, replace=False,
                                             p=available_nodes_to_trace_list_p/sum(available_nodes_to_trace_list_p))
        else:
            sampled_nodes = random.sample(available_nodes_to_trace, nRT)
        for n2t in sampled_nodes:
            if n2t in self.ctd or n2t in self.tlist:
                print('error')
            # try to trace it
            rt = random.random()
            if rt < self.p2bRTs[self.status[n2t]]:
                # Only counts when it is actually traced, as taking up a tracing resource
                if self.tlist_size > 0 and self.tday + self.tlist_size <= self.tmax:
                    self.tlist.add(n2t)
                    self.tlist_remove_time[self.tday + self.tlist_size].append(n2t)
                # it could be it was already in c2t (we allow this), if so remove it
                if n2t in self.c2t:
                    del self.c2t[n2t]
                if self.status[n2t] in self.positive_induced:
                    if n2t not in self.ctd:
                        self.ctd.add(n2t)
                        if n2t not in self.node_trace_depth:
                            self.node_trace_depth[n2t] = 0
                    else:
                        print("And this is impossible")
                    RTpos += 1
                    self.tested_positive_total.add(n2t)
                    # find out contacts
                    self.find_out_contacts(n2t)
                    self.quarantine(n2t)
        self.positives_day.append((CTpos, RTpos))

    def action_on_new_event(self, new_status, n2t):
        # if node in ctd already, it has been traced, contacts discovered and attempted to quarantine, we do not
        # try again
        if self.oracle_type.lower() == 'oracle' and new_status in self.to_trace_instantaneous:
            if n2t not in self.ctd and n2t not in self.c2t and n2t not in self.tlist:
                self.c2t[n2t] = 1
        self.status_in_positive_instantaneous(new_status, n2t)


if __name__ == '__main__':
    test = DefaultTracer(G=None, num_tracers=0, pC2T=1, pQ=1, p2bCTs=1, p2bRTs=1, positive_instantaneous=None,
                         positive_induced=None,
                         pCT=1, tlist_size=0, to_trace_instantaneous=None, oracle_type='aa')
