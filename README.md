## License

Copyright (c) 2024, Marvin Knöller, Jörg Nick

This software is released under GNU General Public License, Version 3.
The full license text is provided in the file LICENSE included in this repository 
or can be obtained from http://www.gnu.org/licenses/

## Content of this project

This is a guide to generate the figures that have been used in the work

**„The temporal domain derivative in inverse acoustic obstacle scattering“**

_By Marvin Knöller and Jörg Nick._

You find all the needed Matlab files in the provided project.

## An overview

To obtain the plots in the numerical examples run IP_RKfor1 for Figure 3, IP_RKfor2 for Figure 4 and IP_RKfor3 for Figure 6.
The visualization of the direct problems is obtained by running visualize_scat1 for Figure 2 and visualize_scat2 for Figure 5.
Figure 1 is simply a sketch, for which we do not provide the code.

# A short explanation of the scripts

All **.m** files have a short description at the beginning of the script. Here, we describe their functionality shortly.

- CalderonCalculusHelmholtz.m : Returns the averaged single layer operator V as a function handle.
- F_op.m : Returns SV^{-1}b, what is required during the convolution quadrature method.
- HelmholtzPotentials.m : Returns the single layer operator S as a function handle.
- ip_RK_for* : Runs the direct and afterward the inverse scattering problem based on the RK-CQ method. 
    Note that the results of IP_RKfor1_double.m are not shown in the work. It is however described in Example 1. 
- PlotIterates.m : Plots exact object and the current iterate given some iteration numbers.
- RKdata.m : Returns the coefficients of the RK method.
- applyRKCQ.m : Applies the RKCQ method. Returns the approximation to the scattering problem u and eventually a usefull density eta.
- center_reg.m & center_regder.m : The second regularization parameter and its Frechet derivative.
- create_RK_sol : Returns the RK solution and eventually a density as well.
- create_triangulation_spline.m : Creates a triangulation around the given spline curve. This is only required for visualization purposes.
- curvepenalty.m & curvepenaltyder.m : The first regularization parameter and its Frechet derivative.
- direct_sound_soft_splineRK.m : Simulates the scattering problem based on splines and RKCQ.
- errorongivenpoints.m : Returns the error.
- generatesplinehandles.m : Generates the periodic, cubic spline.
- get_time_points.m : Returns the time points that are needed for the RKCQ.
- getseed.m & seed.txt : Gets the seed for reproducability of the noisy data simulations.
- minimaldistance.m : Returns the minimal distance of the receivers to the boundary of the obstacle.
- mysignal.m & mysignalp.m : Returns function handles for the function f required to set up the incident wave and its derivative.
- planewaveRK.m & pulseRK.m : Returns the incident wave in the case of a plane wave and a pulse.
- starshapespline.m : Returns the discrete geometry required to describe the obstacles.

The folder Precomputed_results includes all the precomputed results that we obtained in order to generate Figure 3,4 and 6.
    In each folder the simulations 2-5 correspond to noisy data. The outcome is not shown in the work, only described in the last remark.

## Requirements and additional information

The simulations have been carried out on clusters using 32, 64 or 128 cores.
Computations have been carried out using the Matlab 2021a version.
Running one simulation takes 1-2 days. 

The number of time steps and points in space may be decreased, what decreases the compuational time significantly. 
The numerical simulations still yield similar results as in the work.