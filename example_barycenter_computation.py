from test_utils import get_sparsesquare_test, show_sparse_data
from column_generation_for_barycenters import exact_barycenter_colgen_voronoi

##########################################################################
## A) Generate and show the sparse test data

sparse_data = get_sparsesquare_test(k=5,n=5)
show_sparse_data(sparse_data)

##########################################################################
# B) Run our column generation method to obtain the (unweighted) Wasserstein barycenter.
barysol, tuples, tupleweights, optval = exact_barycenter_colgen_voronoi(sparse_data)
show_sparse_data([barysol])
