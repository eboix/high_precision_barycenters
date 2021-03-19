# General multimarginal optimal transport problem solved using column generation.
from pulp import * ## import pulp-or functions
import numpy as np
import time

class MOTProblem:
    def __init__(self, mus, problemname='MOT_Problem', log_file=None):
        """
        mus is a list of the k marginal distributions
        cost_fn takes in a k-tuple and returns the corresponding cost.
        cutting_plane takes in the dual variables, and returns a list of
            violating tuples, if there are any
        problemname is just a string naming the problem
        log_file is the name of the file to which the cost at each iteration and timings are saved.
        """

        self.k = len(mus)
        self.ns = [len(mus[i]) for i in range(self.k)]
        self.cutting_plane_fn = None
        self.cost_fn = None
        self.log_file = log_file

        ## Initialize PuLP problem.
        self.prob = LpProblem(problemname,LpMinimize)    # set up the problem.
        self.obj = LpConstraintVar("obj")   # generate a constraint variable that will be used as the objective
        self.prob.setObjective(self.obj)

        # Add one constraint for each mu_{ij}.
        self.constraintList=[]   # list to save constraint variables in
        for i in range(self.k):
            self.constraintList.append([])
            for j in range(self.ns[i]):
                var=LpConstraintVar("C"+str(i) + "_" + str(j),LpConstraintEQ,mus[i][j])
                self.constraintList[i].append(var)
                self.prob+=var

        self.tups = []
        self.tups_set = set()
        self.TupleVars=[]

    def set_cost_fn(self,cost_fn):
        self.cost_fn = cost_fn
    def set_cutting_plane_fn(self,cutting_plane_fn):
        self.cutting_plane_fn = cutting_plane_fn

    def solveWithCurrentVars(self):
        self.prob.writeLP('prob.lp')
        self.prob.solve()  # start solve

        # print(self.prob.constraints)
        currdualval = []
        curridx = 0
        flatteneddual = [self.prob.constraints[i].pi for i in self.prob.constraints]
        for i in range(self.k):
            currdualval.append([])
            for j in range(self.ns[i]):
                currdualval[i].append(flatteneddual[curridx + j])
            curridx += self.ns[i]
        return currdualval

    def solve(self):

        # Variables for logging purposes
        dual_value_at_each_iteration = []
        cumulative_time_at_each_iteration = []
        starttime = time.time()

        while True:
            duals=self.solveWithCurrentVars()

            # Log current time and current value
            currtime = time.time()
            dual_value_at_each_iteration.append(self.get_objective())
            cumulative_time_at_each_iteration.append(currtime - starttime)


            # Generate new tuples for column generation
            new_tups=self.cutting_plane_fn(duals)
            filtered_new_tups = set(new_tups).difference(self.tups_set);

            if filtered_new_tups:
                for new_tup in filtered_new_tups:
                    self.addTuple(new_tup)
            else:
                if self.log_file:
                    f = open(self.log_file, 'w')
                    for i in range(len(dual_value_at_each_iteration)):
                        f.write('{v:.9f} {t:.9f}\n'.format(v=dual_value_at_each_iteration[i], t=cumulative_time_at_each_iteration[i]));
                    f.close()
                return


    def addTuple(self,tup):
        # Add a new variable in the range [0,1], corresponding to the tuple tup.
        self.tups.append(tup)
        self.tups_set.add(tup)
        tup_cost = self.cost_fn(tup)
        currcolumn = lpSum(self.obj*tup_cost + [self.constraintList[i][tup[i]] for i in range(self.k)] )
        var = LpVariable("Tup"+str(len(self.tups)), 0, 1, LpContinuous, currcolumn)
        self.TupleVars.append(var)

    def get_objective(self):
        return value(self.prob.objective)

    def get_used_tuples(self):
        usedTupList = []
        for i,x in enumerate(self.TupleVars):
            if value(x)>0:
                usedTupList.append((value(x),self.tups[i]))
        return usedTupList
