Files in this directory:

make_synthetic_data.m
- A matlab script that creates synthetic trajectory data.

trajectory_data.txt
- An output of make_synthetic_data.m

bm_skel.m
- The script that runs the nested sampling algorithm. It is currently set to load the data in the file 'trajectory_data.txt'. Note that the 'invprior' fields of the 'models' arrays are currently set to some strange skewed prior, which you probably want to change

bm_results.txt
- A summary file created by a run of bm_skel.m

ns_*.m
- Matlab files that contains the program for running the nested sampling algorithm
