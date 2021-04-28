[![arXiv](https://img.shields.io/badge/arXiv-2006.10415-B31B1B.svg)](https://arxiv.org/abs/2006.10415)
[![MIT Licence](https://badges.frapsoft.com/os/mit/mit.svg?v=103)](https://opensource.org/licenses/mit-license.php)


# Dark Photon Cookbook
Python-3 Code to reproduce the results from our paper arXiv:[2105.XXXXX] "Dark photons: a cookbook"

The code is relatively uncomplicated, just a serious of short python notebooks evaluating the expressions in the paper, and oing some simple Monte Carlo distributions. 
If you are here looking for the code and data needed to make the Dark Photon limit plots, you need to head to [this repo](https://github.com/cajohare/AxionLimits) instead.

If you need any assistance or have any questions contact me at ciaran.aj.ohare@gmail.com

# Requirements
* [`cmocean`](https://matplotlib.org/cmocean/)
* [`numba`](http://numba.pydata.org/)

# Contents
* `data/` - Contains various data files, fluxes, solar models and axion limit data
* `src/` - Main python functions for doing the meat of the analysis
* `notebooks/` - for plotting and doing some extra analysis not found in the main paper
* `plots/` - Plots in pdf or png format

# Examples:
Click to go to the notebook used to make the plot

[<img src="plots/plots_png/LocationDependence.png" width="1000">](https://github.com/cajohare/DarkPhotonCookbook/blob/master/Plot_Bfield.ipynb)

---
