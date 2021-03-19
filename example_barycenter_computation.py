from barycenter_utils import get_sparsesquare_test, show_sparse_data
from barycenter_column_generation import BarycenterProblem

##########################################################################
## A) Generate and show the sparse test data
## Generate k distributions, each uniform on n random points in the unit square.
## sparse_distributions is a list of k distributions. Each distribution is
## represented by a n x 3 numpy array, where each row [x, y, mu] represents a
## support point (x,y) along with its weight mu.
sparse_distributions = get_sparsesquare_test(k=10,n=2)
# show_sparse_data(sparse_distributions) # Draw the distributions

##########################################################################
# B) Run our column generation method to obtain the Wasserstein barycenter.
# method = 'powerdiagram', or 'naive'.
# The powerdiagram method runtime per iteration grows as poly(n,k).
# The naive method runtime grows as n^k.
bary_prob=BarycenterProblem(sparse_distributions,method='powerdiagram')
bary_prob.solve()


print('final_cost of barycenter:',bary_prob.get_objective())
used_tups=bary_prob.get_used_tuples()
print('used_tups for barycenter:',used_tups)
barysol = bary_prob.get_solved_barycenter()

# show_sparse_data([barysol]) # Draw the barycenter solution
