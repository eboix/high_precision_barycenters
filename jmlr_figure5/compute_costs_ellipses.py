import os
import os.path
import numpy as np
import pickle
import sys
sys.path.append('../')
from barycenter_utils import exact_cost_lemon

def compute_barycenter_costs(barycenter_files, input_distributions):

    fs = []
    vals = []

    for f in barycenter_files:
        curr_barycenter = pickle.load(open(f, 'rb'))
        curr_barycenter[:,0] = 1 - curr_barycenter[:,0]

        totcst = 0
        for input_distribution in input_distributions:
            totcst += exact_cost_lemon(input_distribution, curr_barycenter, lemon_solver_location='../lemon_solver/LemonNetworkSimplex')
        vals.append((totcst, f))

    vals.sort()
    return vals
