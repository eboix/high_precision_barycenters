
# Code for Figure 3

This folder contains code to reproduce Figure 3 from the paper [Altschuler, Boix-Adsera JMLR 2021](https://jmlr.org/papers/v22/20-588.html). Figure 3 shows that our algorithm can compute exact barycenters at previously intractable sizes.

![image]('paperfig3.png')

#### How to use
After completing basic installation outlined in [README.md](../README.md), perform the following extended installation.
1. Install MATLAB with the Image Processing Toolbox.

2. Set the variable `matlab_path` in the `barycenter_methods_comparison.py` file to point to your installation of Matlab (e.g., `matlab_path = /Applications/MATLAB_R2020b.app/bin/matlab'`).

3. Install [Lemon](https://lemon.cs.elte.hu/trac/lemon)

4. Compile the fast network simplex solver for pairwise OT that we use to compute the exact cost of a barycenter. This code is a slight modification of [Dong et al.](https://github.com/twistedcubic/fast_ot):
```
cd lemon_solver
g++ -std=c++11 -O3 LemonNetworkSimplex.cpp -o LemonNetworkSimplex -lemon
cd ..
```

Then run `generate_figure_3.py` with Python to compare against
* MAAIPM ([Ge et al. 2019](https://papers.nips.cc/paper/2019/hash/0937fb5864ed06ffb59ae5f9b5ed67a9-Abstract.html))
* Iterated Bregman Projection ([Solomon et al. 2015](https://dl.acm.org/doi/10.1145/2766963))
