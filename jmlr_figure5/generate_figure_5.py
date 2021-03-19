from draw_ellipse_sol import draw_data_set
from compute_costs_ellipses import compute_barycenter_costs
import pickle

# input_distributions is a list of 10 images of concentric ellipses.
# Each image is represented in sparse format as a n x 3 numpy array,
# where each row [x, y, mu] represents a 2-D point (x,y) and probability mass mu on this point.
input_distributions = pickle.load(open('data/input_ellipse_data.pkl', 'rb'))

# This list of files contains the barycenters of the ellipse data that are computed by different methods,
# as outlined in the paper. Each of these files contains a n x 3 numpy array,
# where each row [x, y, mu] represents a 2-D point (x,y) and probability mass mu on this point.
barycenter_data = ['ours_exact.pkl', 'maaipm.pkl', 'debiased_eps002.pkl', 'ibp_eps002.pkl', 'frank_wolfe.pkl']

# Compute costs of all barycenters and draw them in the outputfigs folder.
draw_data_set('data/input_ellipse_data.pkl')
for f in barycenter_data:
    draw_data_set('data/' + f)
    pass
print('Data drawn to outputfigs folder')

try:
    all_costs = compute_barycenter_costs(['data/' + f for f in barycenter_data], input_distributions)

    print('\nCost for each barycenter: ', all_costs)
except Exception as e:
    print('To compute the costs of the barycenters you need to compile the lemon network simplex solver.')
    print(e)


## If you would like to generate our exact ellipse barycenter from scratch, you
# can run the following code. It takes some time to run, as our column generation
# implementation is not optimized.
def solve_ellipse_barycenter_problem():
    import sys
    sys.path.append('../')
    from barycenter_column_generation import BarycenterProblem

    input_distributions = pickle.load(open('data/input_ellipse_data.pkl', 'rb'))
    bary_prob=BarycenterProblem(input_distributions,method='powerdiagram')
    bary_prob.solve()
    barysol = bary_prob.get_solved_barycenter()

    pickle.dump(barysol, open("data/ours_exact_just_computed.pkl", "wb"))

# solve_ellipse_barycenter_problem()

print('\n\nP.S. If you want to test your own algorithm on the ellipses dataset use the "solve_ellipse_barycenter_problem()" method in this file as a template.')
