from MOT_column_generation import MOTProblem
from functools import partial
from barycenter_separation_oracle import get_tuples_to_check_from_power_diagram_intersections, get_tuple_cost_barycenter
import itertools
import numpy as np


class BarycenterProblem(MOTProblem):
    def __init__(self, distributions, method='powerdiagram', problemname='Barycenter_Problem',log_file=None):
        # distributions is a k-length list of numpy n_i x 3 arrays, each row is [x, y, mu]
        k = len(distributions)
        ns = [len(distributions[i]) for i in range(k)]

        mus = []
        self.supportpoints = []
        for i in range(k):
            mus.append(distributions[i][:,2])
            mus[i] = mus[i] / np.sum(mus[i])
            self.supportpoints.append(distributions[i][:,0:2])

        super(BarycenterProblem, self).__init__(mus,problemname,log_file)

        if method == 'powerdiagram':
            # Initialize cost function and cutting plane method for power diagram method.
            # Cost is a function of a k-tuple
            self.cost_fn = partial(get_tuple_cost_barycenter,self.supportpoints)

            # Cutting plane is a function of the dual weights
            self.cutting_plane_fn = partial(get_tuples_to_check_from_power_diagram_intersections,self.supportpoints)

            # Initialize tuples in some arbitrary manner so that a coupling is
            # possible that has the desired marginals and is supported on these
            # tuples
            init_tups = BarycenterProblem.get_waterfilling_initialization_tuples(mus)
            for tup in init_tups:
                self.addTuple(tup)

        elif method == 'naive':
            # Initialize cost function and cutting plane method for naive method.

            # Cost is a function of a k-tuple
            # partial(get_tuple_cost,self.supportpoints)
            self.cost_fn = partial(get_tuple_cost_barycenter,self.supportpoints)

            # Cutting plane is a function of the dual weights
            # The cutting plane function should return an empty list in this naive
            # implementation.
            self.cutting_plane_fn = cutting_plane_naive

            # Since we are running the naive algorithm, initialize so that all of the
            # tuples are there to begin with.
            tupranges = []
            for i in range(k):
                tupranges.append(range(ns[i]))
            for tup in itertools.product(*tupranges):
                self.addTuple(tup)

        else:
            assert(0, 'Error: Method name should be "powerdiagram" or "naive"')



    def get_waterfilling_initialization_tuples(mus):
        print("Waterfilling initialization")
        k = len(mus)
        ns = [len(mus[i]) for i in range(k)]
        curr_tup_idx = []
        for i in range(k):
            curr_tup_idx.append(0)

        waterfillingschedule = []
        for i in range(k):
            currfill = 0
            for j in range(ns[i]):
                currfill += mus[i][j]
                waterfillingschedule.append((currfill,i,j))
        waterfillingschedule.sort()
        # print('waterfillingschedule', waterfillingschedule)

        tup_init_list = []
        currmass = 0
        curr_indices = [0 for i in range(k)]
        tup_init_list.append(tuple(curr_indices))
        for i in range(len(waterfillingschedule)):
            curridx = waterfillingschedule[i][1]
            currreplacement = waterfillingschedule[i][2]
            curr_indices[curridx] = min(currreplacement + 1, ns[curridx] - 1)
            if i < len(waterfillingschedule) - 1:
                if waterfillingschedule[i+1][0] == waterfillingschedule[i][0]:
                    continue
            tup_init_list.append(tuple(curr_indices))

        tup_init_list = list(set(tup_init_list))
        return tup_init_list

    def get_solved_barycenter(self):
        used_tups = self.get_used_tuples()
        supp_size = len(used_tups)
        barysol = np.zeros((supp_size, 3))
        k = len(self.supportpoints)

        barysol = []
        for tup in used_tups:
            tup_mass = tup[0]
            tup_indices = tup[1]
            supp_pt = np.zeros(2)
            for i in range(k):
                supp_pt += self.supportpoints[i][tup_indices[i],:]
            supp_pt = supp_pt / k
            weighted_supp_pt = np.zeros(3)
            weighted_supp_pt[0:2] = supp_pt
            weighted_supp_pt[2] = tup_mass
            barysol.append(weighted_supp_pt)
        print(barysol)
        return np.asarray(barysol)

def cutting_plane_naive(dualwts):
    return []
