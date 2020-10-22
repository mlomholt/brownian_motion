Files in this directory:

make_synthetic_data.m
- A matlab script that creates synthetic trajectory data.

trajectory_data.txt
- An output of make_synthetic_data.m

bm_skel.m
- A Matlab script that runs the nested sampling algorithm. It is currently set to load the data in the file 'trajectory_data.txt'.

bmwmn_skel.m
- Similar to bm_skel.m, except that this version includes models with measurement noise

bmwmn_logl.m
- A matlab function that calculates the log-likelihood for brownian motion with drift and measurement noise

bm_results.txt
- A summary file created by a run of bm_skel.m

bmwmn_results.txt
- A summary file created by a run of bmwmn_skel.m

ns_*.m
- Matlab files that contains the program for running the nested sampling algorithm
