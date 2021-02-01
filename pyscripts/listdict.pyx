# the original code is from https://github.com/springer-math/Mathematics-of-Epidemics-on-Networks
# this was re-written in cython for computational efficiency
import random
from collections import defaultdict, Counter

cdef class _ListDict_(object):
    cdef public dict item_to_position
    cdef public list items
    cdef public bint weighted
    cdef public object weight
    cdef public long double max_weight
    cdef public long double _total_weight
    cdef public long double max_weight_count
    def __init__(self, bint weighted=False):
        self.item_to_position = {}
        self.items = []

        self.weighted = weighted
        if self.weighted:
            self.weight = defaultdict(int)  # presume all weights positive
            self.max_weight = 0
            self._total_weight = 0
            self.max_weight_count = 0

    def __len__(self):
        return len(self.items)

    def __contains__(self, object item):
        return item in self.item_to_position

    cpdef _update_max_weight(self):
        C = Counter(self.weight.values())  # may be a faster way to do this, we only need to count the max.
        self.max_weight = max(C.keys())
        self.max_weight_count = C[self.max_weight]

    cpdef insert(self, object item, long double weight=0):
        r'''
        If not present, then inserts the thing (with weight if appropriate)
        if already there, replaces the weight unless weight is 0

        If weight is 0, then it removes the item and doesn't replace.

        WARNING:
            replaces weight if already present, does not increment weight.


        '''
        if self.__contains__(item):
            self.remove(item)
        if weight != 0:
            self.update(item, weight_increment=weight)

    cpdef update(self, object item, long double weight_increment=-1):
        r'''
        If not present, then inserts the thing (with weight if appropriate)
        if already there, increments weight

        WARNING:
            increments weight if already present, cannot overwrite weight.
        '''
        if weight_increment > 0:  # will break if passing a weight to unweighted case
            if weight_increment > 0 or self.weight[item] != self.max_weight:
                self.weight[item] = self.weight[item] + weight_increment
                self._total_weight += weight_increment
                if self.weight[item] > self.max_weight:
                    self.max_weight_count = 1
                    self.max_weight = self.weight[item]
                elif self.weight[item] == self.max_weight:
                    self.max_weight_count += 1
            else:  # it's a negative increment and was at max
                self.max_weight_count -= 1
                self.weight[item] = self.weight[item] + weight_increment
                self._total_weight += weight_increment
                self.max_weight_count -= 1
                if self.max_weight_count == 0:
                    self._update_max_weight
        elif self.weighted:
            raise Exception('if weighted, must assign weight_increment:', weight_increment)

        if item in self:  # we've already got it, do nothing else
            return
        self.items.append(item)
        self.item_to_position[item] = len(self.items) - 1

    cpdef remove(self, object choice):
        if not self.__contains__(choice):
            return
        position = self.item_to_position.pop(choice)
        last_item = self.items.pop()
        if position != len(self.items):
            self.items[position] = last_item
            self.item_to_position[last_item] = position

        if self.weighted:
            weight = self.weight.pop(choice)
            self._total_weight -= weight
            if weight == self.max_weight:
                # if we find ourselves in this case often
                # it may be better just to let max_weight be the
                # largest weight *ever* encountered, even if all remaining weights are less
                #
                self.max_weight_count -= 1
                if self.max_weight_count == 0 and len(self) > 0:
                    self._update_max_weight()

    cpdef object choose_random(self):
        # r'''chooses a random node.  If there is a weight, it will use rejection
        # sampling to choose a random node until it succeeds'''
        if self.weighted:
            while True:
                choice = random.choice(self.items)
                if random.random() < self.weight[choice] / self.max_weight:
                    break
            # r = random.random()*self.total_weight
            # for item in self.items:
            #     r-= self.weight[item]
            #     if r<0:
            #         break
            return choice

        else:
            return random.choice(self.items)

    cpdef object random_removal(self):
        r'''uses other class methods to choose and then remove a random node'''
        choice = self.choose_random()
        self.remove(choice)
        return choice

    cpdef long double total_weight(self):
        if self.weighted:
            return self._total_weight
        else:
            return len(self.items)

    cpdef update_total_weight(self):
        self._total_weight = 0
        for item in self.items:
            self._total_weight += self.weight[item]
