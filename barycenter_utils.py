
# The following file provides utility methods to generate
# Each sparse data distribution is represented by a list of (x,y,weight) tuples
# Each dense data distribution is represented by an image

import numpy as np
import os
import pylab
import matplotlib.image

###############
# Generate sparse testing instances

def get_sparsedisc_test(k,n):
    """
    Return k sparse n-point distributions on the unit disk
    """
    supports = []
    for i in range(k):
        angle = (2 * np.pi) * np.random.random(n)
        unif_pts = np.stack([np.cos(angle), np.sin(angle)], axis = 1)
        unif_pts *= np.sqrt(np.random.random(n))[:,None]
        pt_wts = np.ones((n,1))
        new_dist = np.concatenate((unif_pts, pt_wts), axis=1)
        supports.append(new_dist)

    return supports

def get_sparsesquare_test(k,n):
    """
    Return k sparse n-point distributions on the unit square
    """
    supports = []
    for i in range(k):
        unif_pts = np.random.random((n,2)) * 2 - 1
        pt_wts = np.ones((n,1))
        new_dist = np.concatenate((unif_pts, pt_wts), axis=1)
        supports.append(new_dist)

    return supports

def get_noisygrid_test(k,n):
    """
        Generate a grid on n points & return k noisy perturbations of it.
        Contained in unit square.
    """
    side_length = math.ceil(math.sqrt(n))
    i = 0
    j = 0
    t = 0
    grid_pts = []
    while t < n:
        grid_pts.append([i, j])
        j += 1
        if j == side_length:
            j = 0
            i += 1
        t += 1
    grid_pts = np.asarray(grid_pts)
    grid_pts = (grid_pts / (side_length-1))*2 - 1

    grid_pts = grid_pts * 0.6

    grid_perturb = np.random.random() * 0.2;
    grid_pts = grid_pts + grid_perturb

    supports = []
    for i in range(k):
        perturb = (2*np.random.random((n,2)) - 1) * 0.2
        pt_wts = np.ones((n,1))
        new_dist = np.concatenate((grid_pts + perturb, pt_wts), axis=1)
        supports.append(new_dist)

    return supports

def get_translatedgrid_test(k,n):
    """
        Generate a grid on n points & return k noisy perturbations of it.
        Contained in unit square.
    """
    side_length = math.ceil(math.sqrt(n))
    i = 0
    j = 0
    t = 0
    grid_pts = []
    while t < n:
        grid_pts.append([i, j])
        j += 1
        if j == side_length:
            j = 0
            i += 1
        t += 1
    grid_pts = np.asarray(grid_pts)
    grid_pts = (grid_pts / (side_length-1))*2 - 1

    grid_pts = grid_pts * 0.7

    grid_perturb = np.random.random() * 0.1;
    grid_pts = grid_pts + grid_perturb

    supports = []
    for i in range(k):
        perturb = (2*np.random.random()-1) * 0.2
        pt_wts = np.ones((n,1))
        new_dist = np.concatenate((grid_pts + perturb, pt_wts), axis=1)
        supports.append(new_dist)

    return supports



###############
# Generate dense testing instances
#

def rebin(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)

def get_dense_img_test(test_set='simple', k=None):

    if test_set == 'simple':

        img_shape = (16,16)
        f1 = 1 - pylab.imread('../data/redcross')[:, :, 2]
        f2 = 1 - pylab.imread('../data/duck')[:, :, 2]
        f3 = 1 - pylab.imread('../data/heart')[:, :, 2]
        f4 = 1 - pylab.imread('../data/tooth')[:, :, 2]

        A = []
        f1 = f1 / np.sum(f1)
        f2 = f2 / np.sum(f2)
        f3 = f3 / np.sum(f3)
        f4 = f4 / np.sum(f4)

        f1 = rebin(f1, img_shape)
        f2 = rebin(f2, img_shape)
        f3 = rebin(f3, img_shape)
        f4 = rebin(f4, img_shape)

        A.append(f1)
        A.append(f2)
        A.append(f3)
        A.append(f4)
        A = np.array(A)

        return A
    elif test_set == 'mnist1':
        img_nums = [3, 6, 8, 14, 23, 24, 40, 59, 67, 70, 72, 77, 78, 99, 102, 105, 112, 113, 124, 128, 134, 152, 174, 177, 184, 200, 201, 205, 208, 211, 224, 231, 248, 251, 269, 270, 276, 290, 309, 310, 315, 345, 351, 355, 357, 358, 366, 382, 394, 397, 398, 406, 408, 416]
        A = []
        if k is None:
            k = 1
        for i in range(k):
            f1 = pylab.imread('../data/1/' + str(img_nums[i]) + '.png')
            f1 = f1 / np.sum(f1)
            f1 = rebin(f1, (14,14))
            A.append(f1)
        return A
    elif test_set == 'mnist2':
        img_nums = [5, 16, 25, 28, 76, 82, 109, 117, 120, 122, 143, 159, 161, 171, 178, 180, 187, 189, 190, 199, 213, 220, 233, 252, 253, 262, 268, 277, 308, 317, 318, 325, 339, 347, 360, 365, 375, 378, 381, 385, 390, 391]
        A = []
        if k is None:
            k = 1
        for i in range(k):
            f1 = pylab.imread('../data/2/' + str(img_nums[i]) + '.png')
            f1 = f1 / np.sum(f1)
            f1 = rebin(f1, (28,28))
            A.append(f1)
        return A
    else:
        raise Exception('Test set ' + str(test_set) + ' not supported.')

#############
# Convert between sparse and dense data representation
#

def get_bounding_box(sparse_data):

    xmin, ymin = math.inf, math.inf
    xmax, ymax = -math.inf, -math.inf
    k = len(sparse_data)
    for i in range(k):
        for j in range(sparse_data[i].shape[0]):
            x, y = sparse_data[i][j,0], sparse_data[i][j,1]
            xmin = min(xmin, x)
            ymin = min(ymin, y)
            xmax = max(xmax, x)
            ymax = max(ymax, y)
    assert(xmin <= xmax)
    assert(ymin <= ymax)
    return ((xmin, ymin), (xmax, ymax))

def get_dense_data_representation(sparse_data, eps, bounding_box=None):

    # Return an image for each element of sparse_data.
    # Each image will be 1/\eps x 1/\eps.


    if bounding_box is None:
        bounding_box = get_bounding_box(sparse_data)

    # The bounding box's parameters.
    xmin = bounding_box[0][0]
    ymin = bounding_box[0][1]
    xmax = bounding_box[1][0]
    ymax = bounding_box[1][1]

    side_length = int(1/eps)
    assert(side_length >= 1)

    k = len(sparse_data)
    imgs = []
    for i in range(k):
        curr_img = np.zeros((side_length, side_length))
        for j in range(sparse_data[i].shape[0]):
            x,y,w = sparse_data[i][j,0], sparse_data[i][j,1], sparse_data[i][j,2]
            r, c = 0,0
            if xmin != xmax:
                c = max(0, min(int((x - xmin) / (xmax - xmin) * side_length), side_length-1))
            if ymin != ymax:
                r = side_length-1-max(0, min(int((y - ymin) / (ymax - ymin) * side_length), side_length-1))
            curr_img[r,c] += w
        imgs.append(curr_img)

    return imgs

def save_dense_data_representation(dense_data, fprefix):
    for i, img in enumerate(dense_data):
        matplotlib.image.imsave(fprefix + '_' + str(i+1) + '.png', 1-img, cmap='binary')

def save_sparse_data_representation(sparse_data, fprefix):

    f = open(fprefix + '.txt', 'w')
    for i, sparse_datum in enumerate(sparse_data):
        for j in range(sparse_datum.shape[0]):
            f.write(str(sparse_datum[j,0]))
            f.write(' ')
            f.write(str(sparse_datum[j,1]))
            f.write(' ')
            f.write(str(sparse_datum[j,2] / np.sum(sparse_datum[:,2])))
            f.write('\n')
    f.close()

def get_sparse_data_representation(dense_data, bounding_box=None):

    # Return a list of points for each element of dense_data.
    # Each point will have a weight corresponding to the weight in the image.

    k = len(dense_data)
    if bounding_box is None:
        bounding_box = ((0,0),(1,1))

    # The bounding box's parameters.
    xmin = bounding_box[0][0]
    ymin = bounding_box[0][1]
    xmax = bounding_box[1][0]
    ymax = bounding_box[1][1]

    pt_lists = []
    for i in range(k):
        curr_pt_list = []
        curr_img = dense_data[i]

        num_rows, num_cols = curr_img.shape
        for r in range(num_rows):
            for c in range(num_cols):
                w = curr_img[r,c]
                x = xmin + (xmax - xmin) * c / num_cols
                y = ymin + (ymax - ymin) * (num_rows - r-1) / num_rows
                assert(w >= 0)
                if w > 0:
                    curr_pt_list.append((x,y,w))
        curr_pts = np.asarray(curr_pt_list)
        pt_lists.append(curr_pts)
    return pt_lists

def show_sparse_data(sparse_data):

    if type(sparse_data) != type(list()):
        sparse_data = [sparse_data]
    for pt_list in sparse_data:
        pylab.scatter(pt_list[:,0], pt_list[:,1], c=pt_list[:,2])
        pylab.show()

def show_dense_data(dense_data):
    for img in dense_data:
        pylab.imshow(img)
        pylab.show()

# Compute exact cost between two sparse images, by solving the LP (Using the LEMON solver).
def exact_cost_lemon(sparse_datum1, sparse_datum2, lemon_solver_location='./lemon_solver/LemonNetworkSimplex'):
    # sparse_datum1 and sparse_datum2 are n1 x 3 and n2 x 3 numpy arrays, respectively
    n1 = sparse_datum1.shape[0]
    n2 = sparse_datum2.shape[0]

    print('Current directory is: ', os.getcwd())
    print('Looking here for the Lemon solver: ', lemon_solver_location)

    if not os.path.isfile(lemon_solver_location):
        assert 0, 'ERROR: Cannot find LemonNetworkSimplex executable. Please compile it (see either jmlr_figure3/README.md or jmlr_figure5/README.md for instructions).'
    print('Lemon solver found!')

    f = open('curr_cost_problem', 'w')
    f.write(str(n1) + ' ' + str(n2))
    f.write('\n')
    norm_const1 = np.sum(sparse_datum1[:,2])
    norm_const2 = np.sum(sparse_datum2[:,2])
    for i in range(n1):
        f.write(str(sparse_datum1[i,0]))
        f.write(' ')
        f.write(str(sparse_datum1[i,1]))
        f.write(' ')
        f.write(str(sparse_datum1[i,2] / norm_const1))
        f.write(' ')
    f.write('\n')
    for j in range(n2):
        f.write(str(sparse_datum2[j,0]))
        f.write(' ')
        f.write(str(sparse_datum2[j,1]))
        f.write(' ')
        f.write(str(sparse_datum2[j,2] / norm_const2))
        f.write(' ')
    f.write('\n')
    f.close()

    os.system(lemon_solver_location + ' < curr_cost_problem > curr_sol_file')
    f = open('curr_sol_file', 'r')
    output = f.read().split()
    f.close()
    os.system('rm curr_cost_problem')
    os.system('rm curr_sol_file')
    val = float(output[0])
    return val
