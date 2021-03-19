import pickle
import numpy as np
import matplotlib.pyplot as plt
import math

def process_img(img,fname=None):
    # Image is an n x 3 numpy array, where each row [x,y,mu] represents
    # a support point at (x,y) of weight mu.
    x = img[:,0]
    y = img[:,1]
    x = 1 - x
    weights = img[:,2]
    rangex = [0,1]
    rangey = [0,1]

    totwt = np.sum(weights)
    weights = weights / totwt

    plt.figure()
    for i in range(len(x)):
        wt_scaling = math.sqrt(5)
        curr_circle = plt.Circle((y[i],x[i]),math.sqrt(weights[i]) * 0.02 * wt_scaling,facecolor='k',edgecolor=None)
        plt.gcf().gca().add_artist(curr_circle)

        plt.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        left=False,
        right=False,
        labelleft=False,
        labelbottom=False) # labels along the bottom edge are off
    plt.axis('scaled')
    plt.xlim(rangex)
    plt.ylim(rangey)
    if fname:
        plt.savefig('outputfigs/' + fname.split('/')[-1] + '.pdf',bbox_inches="tight")
        plt.close()
    else:
        plt.show()

def draw_data_set(filename):
    print('Drawing: ' + filename + ' to outputfigs folder')
    imgs = pickle.load(open(filename, 'rb'))
    if type(imgs) != list:
        imgs = [imgs]
    for i, img in enumerate(imgs):
        process_img(img,fname=filename + '-' + str(i))
