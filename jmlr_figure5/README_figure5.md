# Code for Figure 5

This folder contains code to reproduce Figure 5 from the paper [Altschuler, Boix-Adsera JMLR 2021](https://jmlr.org/papers/v22/20-588.html). Figure 5 shows that the high-precision barycenters computed by our algorithm yield sharper visualizations. This is demonstrated on a standard benchmark dataset of images of concentric ellipses. The raw images to be averaged, barycenter images produced by state-of-the-art algorithms, and the barycenter image produced by our algorithm can be found in this folder.

#### How to use
After completing basic installation outlined in [README.md](../README.md), run `generate_figure_5.py` with Python:
1. Draws the input data and barycenters computed by different methods to the `outputfigs` folder.

2. (OPTIONAL) Calculates the cost of the barycenter computed by each method. This second step requires the network simplex in `../lemon_solver` to be compiled:
   * Install [Lemon](https://lemon.cs.elte.hu/trac/lemon)

   * Compile the fast network simplex solver for pairwise OT that we use to compute the exact cost of a barycenter. This code is a slight modification of [Dong et al.](https://github.com/twistedcubic/fast_ot):
  ```
  cd ../lemon_solver
  g++ -std=c++11 -O3 LemonNetworkSimplex.cpp -o LemonNetworkSimplex -lemon
  ```

The program `generate_figure_5.py` also contains instructions to help you compare your own barycenter algorithm on the ellipses dataset.

#### Contents of `data` folder
The `data` folder contains:
* `ellipses_data.pkl`: the input dataset (nested ellipses)
* `ours_exact.pkl`: the exact barycenter computed by our method ([Altschuler, Boix-Adsera 2021](https://jmlr.org/papers/v22/20-588.html))

As well as the approximate barycenters computed by comparison methods:
* `maaipm.pkl`: the barycenter computed by MAAIPM ([Ge et al. 2019](https://papers.nips.cc/paper/2019/hash/0937fb5864ed06ffb59ae5f9b5ed67a9-Abstract.html))
* `debiased_eps002.pkl`: the barycenter computed by Debiased Sinkhorn ([Janati et al. 2020](http://proceedings.mlr.press/v119/janati20a.html))
* `ibp_eps002.pkl`: the barycenter computed by Iterated Bregman Projection ([Solomon et al. 2015](https://dl.acm.org/doi/10.1145/2766963))
* `frank_wolfe.pkl`: the barycenter computed by Frank-Wolfe ([Luise et al. 2019](https://papers.nips.cc/paper/2019/hash/9f96f36b7aae3b1ff847c26ac94c604e-Abstract.html))
