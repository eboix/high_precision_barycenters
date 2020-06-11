from barycenter_separation_oracle import get_tuples_to_check_from_power_diagram_intersections
import numpy as np
from cylp.cy import CyClpSimplex
from cylp.py.modeling.CyLPModel import CyLPArray
import time
from cylp.py.pivots import PositiveEdgePivot
from cylp.py.pivots.WolfePivot import WolfePivot
from cylp.py.pivots import LIFOPivot
from cylp.py.pivots import MostFrequentPivot
from cylp.py.pivots import DantzigPivot

import matplotlib.pyplot as plt
import math

def get_mot_primal_solution(supports, basic_inds):
    k = len(supports)
    s = CyClpSimplex()
    num_inds = len(basic_inds)

    x = s.addVariable('x', num_inds)
    s += 0 <= x 
    for i in range(k):
        num_pts_in_curr = len(supports[i])
        a = []
        for j in range(num_pts_in_curr):
            a.append([])
        for idx, ind in enumerate(basic_inds):
            a[ind[i]].append(idx)

        for j in range(num_pts_in_curr):
            s += x[np.asarray(a[j], dtype='int32')].sum() == supports[i][j,2]


    costs = np.zeros(num_inds)
    centroids = np.zeros((num_inds,2))
    for i in range(k):
        centroids += supports[i][basic_inds[:,i],0:2]
    centroids = centroids / k
    for i in range(k):
        costs += np.sum(np.square(supports[i][basic_inds[:,i], 0:2] - centroids),axis=1)

    obj = CyLPArray(costs)
    s.objective = obj * x
    s.primal()
    primsol = s.primalVariableSolution['x']
    usedcentroids = centroids[primsol != 0,:]
    usedtuples = basic_inds[primsol != 0,:]
    usedweights = np.reshape(primsol[primsol != 0], (usedcentroids.shape[0], 1))
    barysol = np.concatenate((usedcentroids, usedweights), axis=1)

    optval = math.inf
    if s.getStatusCode() in [0,3]:
        optval = s.objectiveValue

    return barysol, usedtuples, usedweights, optval


def exact_barycenter_sparse_colgen(supports, tuple_generation_function, log_file=None):
    # supports is a k-length list of numpy n_i x 3 arrays, each row is [x, y, mu]
    # tuple_generation function is a cutting plane oracle for the exact barycenter problem
    # log_file is an optional output file to write the error vs. time at each step
    # 
    # returns (barysol, usedtuples, usedweights, optval), where
    # barysol is the barycenter in sparse format (numpy N x 3 array, where each row is [x, y, mu])
    # usedtuples is a N x k array of tuples of points used to construct the support of the barycenter
    # usedweights is a N x 1 array of weights associated with each tuple in usedtuples
    # optval is the cost of the barycenter

    starttime = time.time()
    k = len(supports)

    # Preprocessing
    M = 0
    num_pts_in_diagram = []
    num_pts = 0
    pt_locs = []
    mus = []
    for i in range(k):
        curr_num_pts = supports[i].shape[0]
        num_pts_in_diagram.append(num_pts)
        num_pts += curr_num_pts

        pt_locs.append(supports[i][:,0:2])

        normalized_wts = supports[i][:,2]
        normalized_wts = normalized_wts / np.sum(normalized_wts)
        mus.extend(normalized_wts.tolist())

        M = max(M, np.max(np.abs(supports[i][:,0:2])))
        M = max(M, np.max(normalized_wts))

    num_pts_in_diagram.append(num_pts)
    mus = np.asarray(mus)

    # Set up the LP (dual MOT)
    s = CyClpSimplex()
    p = s.addVariable('p', num_pts)
    s += p <= 3*M

    curr_tuple_constraints_set = set()

    obj = CyLPArray(mus)
    s.objective = - obj * p

    dual_value_at_each_iteration = []
    cumulative_time_at_each_iteration = []

    # Run the cutting plane method on the dual MOT problem
    while True:
        s.dual()
        
        currtime = time.time()
        dual_value_at_each_iteration.append(s.objectiveValue)
        cumulative_time_at_each_iteration.append(currtime - starttime)

        prim_sol = s.primalVariableSolution['p']
        wts = [prim_sol[num_pts_in_diagram[i]:num_pts_in_diagram[i+1]] for i in range(k)]

        # Find the new tuples (cutting planes) to add, given the current weights.
        tuples_to_check_set = tuple_generation_function(pt_locs, wts)

        new_tuples_to_check_set = tuples_to_check_set.difference(curr_tuple_constraints_set)
        new_tuples_to_check = np.asarray(list(new_tuples_to_check_set))
        num_new_tuples_to_check = len(new_tuples_to_check)

        # If no cutting planes found, we are optimal!
        if num_new_tuples_to_check == 0:
            break

        # Compute the costs associated with each cutting plane at the current point.
        centroids = np.zeros((num_new_tuples_to_check,2))
        for i in range(k):
            centroids += supports[i][new_tuples_to_check[:,i],0:2]
        centroids = centroids / k
        costs = np.zeros(num_new_tuples_to_check)

        for i in range(k):
            costs += np.sum(np.square(supports[i][new_tuples_to_check[:,i], 0:2] - centroids),axis=1)

        curr_tup_wts = np.zeros(num_new_tuples_to_check)
        for i in range(k):
            curr_tup_wts += wts[i][new_tuples_to_check[:,i]]

        # Add a new constraint for each violated cutting plane.
        violations = curr_tup_wts > costs
        viol_costs = costs[violations]
        num_viol_tuples = viol_costs.shape[0]
        viol_pindices = np.zeros((num_viol_tuples, k), dtype=np.int32)
        for i in range(k):
            viol_pindices[:,i] = num_pts_in_diagram[i] + new_tuples_to_check[violations,i]

        for j in range(num_viol_tuples):
            s += p[viol_pindices[j,:]].sum() <= viol_costs[j]

        if num_viol_tuples == 0:
            break

        new_tuples_added = [tuple(x) for x in new_tuples_to_check[violations,:].tolist()]
        curr_tuple_constraints_set.update(new_tuples_added)

        print('Added ', num_viol_tuples, 'new constraints')

#    print('Number of primal variables', len(s.primalVariableSolution['p']))

    # Once the MOT dual has been solved to optimality, extract an MOT primal solution from the active constraints.
    primvarsol = s.primalVariableSolution['p']
    num_nonzero = 0
    num_basic = 0
    basic_inds = []
    for constraint in s.constraints:
        for _, vc in constraint.varCoefs.items():
            inds = vc.nonzero()[1]
            val = 0
            for ind in inds:
                val += primvarsol[ind]
        assert(val <= constraint.upper + 1e-7)
        if abs(val-constraint.upper) < 1e-7:
            num_basic += 1
            basic_inds.append(inds)
    basic_inds = np.asarray(basic_inds).tolist()
    basic_inds = [[x[i] - num_pts_in_diagram[i] for i in range(len(x))] for x in basic_inds]
    basic_inds = np.asarray(basic_inds)
    print('Number of active constraints', num_basic)

    primsol = get_mot_primal_solution(supports, basic_inds)

    # Save the value/wall-time information if a log-file was specified.
    if log_file:
        f = open(log_file, 'w')
        for i in range(len(dual_value_at_each_iteration)):
            f.write('{v:.9f} {t:.9f}\n'.format(v=-dual_value_at_each_iteration[i]/num_pts_in_diagram[-1], t=cumulative_time_at_each_iteration[i]));
        f.close()

    return primsol

def exact_barycenter_colgen_voronoi(supports, log_file=None):
    return exact_barycenter_sparse_colgen(supports, get_tuples_to_check_from_power_diagram_intersections, log_file=log_file)


