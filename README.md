# High-precision Wasserstein Barycenters

This repository contains the code from [Altschuler, Boix-Adsera](https://arxiv.org/abs/2006.08012) for exactly computing Wasserstein barycenters between discrete distributions.

NB: Our focus in writing this code was clarity over speed. Much further optimization can certainly be done.

## Overview

### Main files

* `example_barycenter_computation.py` provides example usage of our code.

* `column_generation_for_barycenters.py` implements the main algorithm for 

* `barycenter_separation_oracle.py` implements the separation oracle for the MOT dual problem.

* `power_diagram_construction.py` contains code to construct power diagrams and their intersection.

* `test_utils.py` contains miscellaneous code to manipulate data.

* `barycenter_methods_comparison.py` contains code to compare to MAAIPM and IBP methods.

### Other files

In the `test_files` directory, the `MAAIPM_code` directory contains the [MAAIPM method code of Ge, Wang, Xiong, and Ye](https://gitlab.com/ZXiong/wasserstein-barycenter), and the `IBP_code` directory contains the [IBP method code of Solomon, de Goes, Peyre, Cuturi, Butscher, Nguyen, Du, Guibas](https://github.com/gpeyre/2015-SIGGRAPH-convolutional-ot). The `lemon_solver` directory contains an exact optimal transport solver implemented by [Dong, Gao, Peng, Razenshteyn, Sawlani](https://github.com/twistedcubic/fast_ot). These files are used by `barycenter_methods_comparison.py` to compare our algorithm to MAAIPM and IBP. The `skgeom_hacked` folder contains a slightly-augmented version of the [skgeom library](https://github.com/scikit-geometry/scikit-geometry), which implements a Python wrapper for the CGAL C++ library. This library is used for manipulating line-segment arrangements.

## Installation (basic, recommended)

This installation allows you to run our barycenter-computation code in your applications (i.e., the code in the `column_generation_for_barycenters.py` folder). In order to run our comparisons against MAAIPM and IBP, the "extended installation" is required.

Requires:
* [CGAL](https://doc.cgal.org/latest/Manual/installation.html) version 5.0.1 or higher

Using [Anaconda](https://docs.anaconda.com/anaconda/install/), create a new environment and install basic dependencies as follows:

```
conda create --name barycenters_environment
conda activate barycenters_environment
conda install -c conda-forge cgal
conda install -c conda-forge matplotlib
conda install -c conda-forge pybind11
conda install -c conda-forge scikit-learn
conda install -c conda-forge opencv
conda install -c conda-forge cvxopt
```

Then build the modified_skgeom library. This is just the [skgeom library](https://github.com/scikit-geometry/scikit-geometry), which wraps CGAL into Python, plus an extra wrapper method that allows for faster construction of line segment arrangements. This installation step is slow. I have submitted a request to add this functionality to the main skgeom library, so hopefully in the future it will be possible to directly install the skgeom library and installation will be simpler and faster.
```
cd skgeom_hacked/scikit-geometry-master
python setup.py build
python setup.py install
cd ../..
```

Then install CyLP, as detailed on the [CyLP website](https://github.com/coin-or/CyLP#cylp). E.g., after installing the binaries for Cbc using coinbrew, run
```
export LD_LIBRARY_PATH=/path/to/coinbrew_folder/build/lib
export COIN_INSTALL_DIR=/path/to/coinbrew_folder/build
pip install cylp
```

## Installation (extended)

This installation allows you to run our comparison tests against MAAIPM and IBP (i.e., the code in `barycenter_methods_comparison.py`). First perform the basic installation above, and then complete it as follows.

Requires:
* MATLAB (to run MAAIPM, IBP implementations)
* [Lemon](https://lemon.cs.elte.hu/trac/lemon)

Compile the fast network simplex solver for pairwise OT that we use to compute the exact cost of a barycenter returned by IBP. This code is a slight modification of the code of [Dong, Gao, Peng, Razenshteyn, Sawlani](https://github.com/twistedcubic/fast_ot):
```
cd lemon_solver
g++ -std=c++11 -O3 LemonNetworkSimplex.cpp -o LemonNetworkSimplex -lemon
cd ..
```

## License

This project is licensed under the LGPL-3 license. See the [LICENSE.md](LICENSE.md) file for details.

