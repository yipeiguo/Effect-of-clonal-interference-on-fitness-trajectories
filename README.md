# Effect-of-clonal-interference-on-fitness-trajectories

This respository contains codes for simulating the evolution of a population under the moran model and for analyzing the simulation data.

System Requirements/Installation:
These codes were ran with MATLAB R2020b.

Instructions for use:

The file 'Script_MoranEvolution_forcluster.m' is the script for running the Moran simulations (it uses the function 'EvolveMoran_Gillispie_v2.m'). Running this script will save the output data in a matlab data file, which is then compared with theoretical predictions using the script `Script_AnalyzeClusterData_v2.m'.

The script 'Script_plotAcrossParams.m' combines data (the output of 'Script_AnalyzeClusterData_v2.m') for all parameters of N and mu for a specified landscape.
The output of this is then used by the following scripts:
- 'Script_fitdata.m' for fitting the trajectories to the equation governing the dynamics of the fitness trajectory in order to extract properties of the fitness landscape
- 'Script_varytscale_theoreticalTrajs.m' for exploring the effect of varying the way the initial fitness gradient is estimated.
