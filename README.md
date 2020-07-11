# Matrix-Calibration-Based Cascaded Channel Estimation for Reconfigurable Intelligent Surface Assisted Multiuser MIMO

This is the simulation code package for the following paper:

Hang Liu, Xiaojun Yuan, and Ying Jun Zhang. "Matrix-calibration-based cascaded channel estimation
for reconfigurable intelligent surface assisted multiuser MIMO," IEEE Journal on Selected Areas in Communications, 2020, doi: 10.1109/JSAC.2020.3007057.

The package, written on Matlab 2014b, reproduces the numerical results in the paper.

## Abstract of Article:

> Reconfigurable intelligent surface (RIS) is envisioned to be an essential component of the paradigm for beyond 5G networks as it can potentially provide similar or higher  array gains with much lower hardware cost and energy consumption compared with the massive multiple-input multiple-output (MIMO) technology. In this paper, we focus on one of the fundamental challenges, namely the channel acquisition, in a RIS-assisted multiuser MIMO system. The state-of-the-art channel acquisition approach in such a system with fully passive RIS elements estimates the cascaded transmitter-to-RIS and RIS-to-receiver channels by adopting excessively long training sequences. To estimate the cascaded channels with an affordable training overhead, we formulate the channel estimation problem in the RIS-assisted multiuser MIMO system as a matrix-calibration based matrix factorization task.  By exploiting the information on the slow-varying channel components and the hidden channel sparsity, we propose a novel message-passing based algorithm to factorize the cascaded channels.  Furthermore, we present an analytical framework to characterize the theoretical performance bound of the proposed estimator in the large-system limit. Finally, we conduct simulations to verify the high accuracy and efficiency of the proposed algorithm.

## How to Use
This package is written on MATLAB 2014b. It includes the following scripts (Also see each file for further documentation):

* __main_Section_VI_A.m__:
This script produces the data for the purple dashed curves in Fig. 5 for given noise power (tau_N); The results are MSE_G_simulation and MSE_S_simulation, which is also stored in DATA/VIA_Simulation.mat. To save the running time, one can set the number of Monte Carlo trials to a small number. To fully recover the plots in Fig. 5, one should set libopt.trials (Line 22) to 5000.

* __main_replica.m__:
This script can produce the solid pruple curve in Fig. 5. Specifically, it provides an iterative algorithm to compute the asymptotic MSEs by computing the fixed-point of eq. (37). The result (MSE_G_ana and MSE_S_ana) are stored in DATA/VIA_Analytical.mat.

* __plot_Fig5.m__: Plot figure as in Fig. 5 based on .mat data files in DATA/.


* __Model_Generation_Library/__: Scripts for generate system models
  * __A_GEN.m__: Generate a BS sampling basis with a unfirom sampling grid;
  * __F_GEN.m__: Generate RIS basis with grids specified by the inputs;
  * __S_GEN.m__: Generate a Bernouli complex Gaussian matrix with unit non-zero variance and specified sparsity;
  * __MSE_Compute.m__: Compute the MSE between the ground-truth and its estimate;
  * __NMSE_CAL.m/NMSE_CAL2.m__: Compute NMSEs by eq. (40).
  
* __DATA/__: Store the numerical results with .mat format.
  
* __MP_Library/__: Scripts for message passing updating
  * __MessagePassing.m/MessagePassing_iteration.m__: Implementation of Algorithm 1;
  * __MPOpt.m__: Class that contains the algorithm control parameters;
  * __EstimIn.m/EstimOut.m/SparseScaEstim.m/CAwgnEstimIn.m/CAwgnEstimOut.m__: Classes contain prior functions for input and output.

* __Replica_Library/replica_iteration.m__: Iterative algorithm for asymptotic MSE computation.


* Codes for simulations in Section VI-B: comming soon
## Referencing

If you in any way use this code for research that results in publications, please cite our original article listed above.

## Acknowledgements

Part of this package is written based on the [GAMPMATLAB package](https://sourceforge.net/projects/gampmatlab/).
