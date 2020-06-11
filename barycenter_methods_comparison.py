import numpy as np
import pylab
import math
from column_generation_for_barycenters import exact_barycenter_colgen_voronoi
import random
import cvxopt
import os
import time
from test_utils import get_sparsedisc_test, get_sparsesquare_test, get_noisygrid_test, get_translatedgrid_test, show_sparse_data, exact_cost_lemon, get_dense_data_representation, save_dense_data_representation, save_sparse_data_representation, get_sparse_data_representation

# Code to run tests of column_generation_for_barycenters against the IBP and MAAIPM methods with fixed support.

def read_val_time_file(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()
    splitlines = [line.split() for line in lines]
    vals = [float(x[0]) for x in splitlines]
    times = [float(x[1]) for x in splitlines]
    return vals, times

def main():

    ## THE FOLLOWING FLAGS TO TRUE/FALSE CONTROL WHICH TESTS ARE DONE
    do_ibp_tests = True
    do_ipm_tests = True
    do_colgen_tests = True
    do_comparison_plots = True
    
    for k in [5]:

        for n in [5]:

            np.random.seed(19238412)
            random.seed(12834981)
            test_type = 'sparsesquare'

            data_identifier = test_type + '_n' + str(n) + 'k' + str(k)

            ##########################################################################
            ## A) Generate the sparse test data
            if test_type == 'sparsesquare':
                sparse_data = get_sparsesquare_test(k=k,n=n)
            elif test_type == 'sparsedisc':
                sparse_data = get_sparsedisc_test(k=k,n=n)
            elif test_type == 'noisygrid':
                sparse_data = get_noisygrid_test(k=k,n=n)
            elif test_type == 'translatedgrid':
                sparse_data = get_translatedgrid_test(k=k,n=n)
            else:
                assert(0)
#            show_sparse_data(sparse_data)

            # And save it for use by IBP and IPM
            data_filename = 'test_files/experiment_data/' + data_identifier

            # Save it as gridded images (for IBP):
            eps_vals = [0.1]
            for i, eps_val in enumerate(eps_vals):
                dense_data = get_dense_data_representation(sparse_data, eps=eps_val, bounding_box=((-1,-1),(1,1)))
                eps_desc = str(eps_val).split('.')[1]
                dense_data_filename = data_filename + 'eps' + eps_desc
                save_dense_data_representation(dense_data, dense_data_filename)

            # Save it as a list of points (for IPM):
            save_sparse_data_representation(sparse_data, data_filename)

            ##########################################################################
            # B) Run our column generation method, collecting an error plot.
            if do_colgen_tests:
                filename_colgen = 'test_files/experiment_results/colgen_results/' + data_identifier + '.txt';
                barysol, tuples, tupleweights, optval = exact_barycenter_colgen_voronoi(sparse_data, log_file=filename_colgen)

            ########################################################################
            ## C) Run IBP method, gridding the space.

            ## 1) Run the IBP code (MATLAB) computing regularized barycenters (different regularizations and epsilons).
            if do_ibp_tests:
                eps_vals = [0.1]
                filter_sizes = [1] # the regularization parameter is eta = (1/eps)^2 / (2 * filter_size^2)
                for i, eps_val in enumerate(eps_vals):
                    eps_desc = str(eps_val).split('.')[1]
                    for filter_size in filter_sizes:
                        ibp_res_prefix = 'test_files/experiment_results/ibp_results/' + data_identifier + 'eps' + eps_desc + 'f' + str(filter_size);
                        ibp_res_file = ibp_res_prefix + '.txt'
                        if os.path.exists(ibp_res_file):
                            print('ALREADY EXISTS: ' + ibp_res_file)
                            continue
                        time.sleep(0.1)
                        matlab_command = 'matlab -nodisplay -r "cd(\'test_files/IBP_code\'); test_barycenters_ibp(' + str(n) + ', ' + str(k) + ',\'' + eps_desc + '\',\'' + test_type + '\',' + str(filter_size) + ',{imgortime:d});exit"';

                        ## i) Assign scores to the IBP barycenters w.r.t. actual point cloud.
                        ## Use the Lemon solver to obtain the exact optimal transport distance between pairs of points.
                        costs_filename = ibp_res_prefix + '_costs.txt'
                        costs_lemon = []
                        if not os.path.exists(costs_filename):
                            os.system(matlab_command.format(imgortime=1)) # Write the barycenter at each iteration to the disk.
                            t = 1
                            while True:
                                curr_filename = ibp_res_prefix + '/b_' + str(t) + '.png'
                                if not os.path.exists(curr_filename):
                                    break
                                if t % 40 == 1 or (t < 80 and t % 10 == 1):
                                    f1 = pylab.imread(curr_filename)
                                    f1 = f1 / np.sum(f1)
                                    sparsedatum_f1 = get_sparse_data_representation([f1], bounding_box=((-1,-1),(1,1)))[0]

                                    print(t)
                                    curr_lemon = 0
                                    for i in range(k):        
                                        print(t,i)
                                        curr_lemon += exact_cost_lemon(sparse_data[i], sparsedatum_f1) / k
                                    costs_lemon.append(curr_lemon)
                                else:
                                    costs_lemon.append(-1)
                                t += 1

                            costs_file = open(costs_filename, 'w')
                            for i in range(len(costs_lemon)):
                                costs_file.write('{v:.9f}\n'.format(v=costs_lemon[i]))
                            costs_file.close()
                        else:
                            print('USING PREVIOUSLY-COMPUTED COSTS FILE: ' + costs_filename)
                            costs_file = open(costs_filename, 'r')
                            costs_lemon = costs_file.readlines()
                            costs_lemon = [float(i) for i in costs_lemon]
                            costs_file.close()
                            
                        # ii) Separately re-run the IBP code to get the timing information.
                        timings_filename = ibp_res_prefix + "_timing.txt";
                        if not os.path.exists(timings_filename):
                            os.system(matlab_command.format(imgortime=2)) # Write the cumulative time to each iteration to the disk.
                        else:
                            print('USING PREVIOUSLY-COMPUTED TIMINGS FILE: ' + timings_filename)
                        timings_file = open(timings_filename, 'r')
                        times = timings_file.readlines()
                        times = [float(i) for i in times]
                        timings_file.close()

                        res_file = open(ibp_res_file, 'w')
                        assert(len(times) == len(costs_lemon))
                        for i in range(len(times)):
                            res_file.write('{v:.9f} {t:.9f}\n'.format(v=costs_lemon[i], t=times[i]));
                        res_file.close()

            ##########################################################################
            ## D) Run MAAIPM method with fixed-support assumption.

            if do_ipm_tests:
                for eps_inv in [10]: # eps_inv x eps_inv-size grid
                    maaipm_res_file = 'test_files/experiment_results/ipm_results/' + data_identifier + 'epsinv' + str(eps_inv) + '.txt'
                    if os.path.exists(maaipm_res_file):
                        print('ALREADY EXISTS: ' + maaipm_res_file)
                        continue
                    time.sleep(0.1)
                    matlab_command = 'matlab -nodisplay -r "cd(\'test_files/MAAIPM_code\'); test_barycenters_maaipm_grid_support(' + str(n) + ', ' + str(k) + ',' + str(eps_inv) + ',\'' + test_type + '\');exit"';
                    
                    os.system(matlab_command)
                    os.system('mv test_files/experiment_results/ipm_results/latest_maaipm_results.txt ' + maaipm_res_file)



    ##############################################################################
    ## E) Create the comparison plots between the algorithms
    if do_comparison_plots:
        print('Comparison plots here')
        markers = ['x', 'o', 'v', '^', '.']*100
        best_val = 10;
        last_time = 0;
        title_text = '';
        legend_names = None;

#        series = ['colgen_results/sparsesquare_n25k10.txt', 'ibp_results/sparsesquare_n25k10eps1f1.txt', 'ibp_results/sparsesquare_n25k10eps1f2.txt', 'ibp_results/sparsesquare_n25k10eps05f1.txt', 'ibp_results/sparsesquare_n25k10eps05f2.txt']
#        series = ['colgen_results/sparsesquare_n25k10.txt', 'ipm_results/sparsesquare_n25k10epsinv10.txt', 'ipm_results/sparsesquare_n25k10epsinv30.txt', 'ipm_results/sparsesquare_n25k10epsinv70.txt']

##        # IPM plot: n = 20, k = 20, noisygrid.
##        series = ['colgen_results/noisygrid_n20k10.txt', 'ipm_results/noisygrid_n20k10epsinv10.txt', 'ipm_results/noisygrid_n20k10epsinv30.txt', 'ipm_results/noisygrid_n20k10epsinv50.txt']

##        # IBP plot: n = 20, k = 20, sparsesquare.
##        series = ['colgen_results/sparsesquare_n20k10.txt', 'ibp_results/sparsesquare_n20k10eps04f1.06.txt', 'ibp_results/sparsesquare_n20k10eps01f7.txt', 'ibp_results/sparsesquare_n20k10eps004f12.txt']
##        legend_names = ['Proposed algorithm', 'IBP 25x25 grid, $\eta=100$', 'IBP 100x100 grid, $\eta=100$', 'IBP 250x250 grid, $\eta=200$']
##        title_text = 'Comparison with IBP'

#        # IPM plot: n = 20, k = 20, sparsesquare.
#        series = ['colgen_results/sparsesquare_n20k10.txt', 'ipm_results/sparsesquare_n20k10epsinv10.txt', 'ipm_results/sparsesquare_n20k10epsinv40.txt', 'ipm_results/sparsesquare_n20k10epsinv70.txt']
#        legend_names = ['Proposed algorithm', 'MAAIPM 10x10 grid', 'MAAIPM 40x40 grid', 'MAAIPM 70x70 grid']
#        title_text = 'Comparison with MAAIPM'

        # Comparison plot: n = 5, k = 5, sparsesquare.
        series = ['colgen_results/sparsesquare_n5k5.txt', 'ipm_results/sparsesquare_n5k5epsinv10.txt', 'ibp_results/sparsesquare_n5k5eps1f1.txt']
        legend_names = ['Proposed algorithm', 'MAAIPM 10x10 grid', 'IBP 10 x 10 grid, filter-size 1']
        title_text = 'Example omparison with MAAIPM and IBP'

        if legend_names is None:
            legend_names = series

        val_data = []
        time_data = []
        
        # Sanitize data (remove places where cost was not computed.)
        for i in range(len(series)):
            temp_vals, temp_times = read_val_time_file('test_files/experiment_results/' + series[i])
            vals = []
            times = []
            for j in range(len(temp_vals)):
                if temp_vals[j] >= 0:
                    vals.append(temp_vals[j])
                    times.append(temp_times[j])
            val_data.append(vals)
            time_data.append(times)

        for i, serie in enumerate(series):
            vals = val_data[i]
            times = time_data[i]
            best_val = min(best_val, np.min(vals)) # Provided that column generation is one of the methods run, best_val will be the true optimum for the problem.
            last_time = max(last_time, np.max(times))

        pylab.figure(figsize=(12,8))
        for i, serie in enumerate(series):
            vals = val_data[i]
            times = time_data[i]
            vals.append(vals[-1])
            times.append(last_time)
            pylab.plot(times, vals-best_val + 1e-18, marker=markers[i], label=legend_names[i], linewidth=3.0)
        pylab.ylim((1e-10 - 1e-11, 10.0))
        pylab.xlim((0.009,10))
        pylab.yscale('log')
        pylab.xscale('log')

        font = {'family' : 'sans-serif',
                'size'   : 30} 
        pylab.rc('font', **font)
        pylab.rcParams['text.usetex'] = True 

        ax = pylab.gca();
        ax.set_ylabel('Suboptimality gap')
        ax.set_xlabel('Time (seconds)')
        ax.set_title(title_text)
        for item in ([ax.xaxis.label, ax.yaxis.label] +
                     ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(30)
            item.set_usetex(True)
        ax.title.set_fontsize(40)
        ax.title.set_usetex(True)

        pylab.legend()
        pylab.show()

if __name__ == '__main__':
    main()

