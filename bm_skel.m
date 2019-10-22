%function bm_skel()  %Use this line if you do want the script to interfere with your workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This skeleton routine carries out Bayesian inference on a trajectory
% A number of parameters may be specified by the user:
%
%  - options.nwalkers is the number of initial samples used in the
%    nested sampling algorithm to sample parameter space.
%
%  - options.stoprat is the ratio between the last added evidence and the
%    total samples evidence at which nested sampling will terminate.
%
%  - options.nsteps is the number steps attempted with the nested sampling MCMC
%    technique to generate a uniformly distributed sample from another sample.
%
%  - the "models" structure specify the functions used for likelihood calculation,
%    sample generation, which are specific 
%    to the system of interest. To use the nested sampling framework on another system, 
%    similar functions must exist and be specified in this structure.
%
% The routine outputs a .txt file with evidence estimations and inferred parameter values.
% The routine writing this file takes the following inputs:
%  - misc.percentiles are the percentiles used for the characterization of the
%    posterier 
%  - misc.labels are the labels assigned to each inferred parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% Load the trajectory data and convert from positions to steps
data=load('trajectory_data.txt');
obs=diff(data);

% There are two models:
% Model 1 is Brownian motion
% Model 2 is Brownian motion with a drift

%Specify the log-likelihood functions for the models
log_normal=@(obs,mu,var) -log(sqrt(2*pi*var))*numel(obs)-sum(sum((obs-mu).^2))/(2*var);
models(1).logl=@(obs,theta) log_normal(obs,0,theta(1)^2);
models(2).logl=@(obs,theta) log_normal(obs,theta(1,2:3),theta(1)^2);

%NOTE: the priors below are probably not realistic, since they use:
inv_skewed=@(u) 2*sqrt(u)-1;   %This inverse cumulative has density (x+1)/2 with -1<x<1

%Specify priors for the models, i.e., the inverse of the cumulative prior distributions
%
%Prior for Model 1
width_sigma=10;
inv_sigma=@(u) width_sigma*(inv_skewed(u)+1)/2;
models(1).invprior=@(u) inv_sigma(u);
%
%Prior for Model 2 (using Model 1's prior for sigma)
width_mu=10;
inv_mu=@(u) width_mu*inv_skewed(u);
models(2).invprior=@(u) [inv_sigma(u(1)), inv_mu(u(2)), inv_mu(u(3))];

%Specify the random u-generators (that provide inputs for the invprior-functions)
models(1).genu=@() rand(1,1);
models(2).genu=@() rand(1,3);

%Specify options
options.Nparfor=1;       % Do not use any parallel computing
models(1).options=options;
models(2).options=options;

%Specify the index for the labels
models(1).labels=[1];
models(2).labels=[1 2 3];

%Labels for the parameters and percentile line in the summary text-file
misc.labels=...
['sigma: ';...
 'mu_x:  ';...
 'mu_y:  '];

%Specify output filename beginnings
misc.data_id = 'bm';

%Tell ns_print to write a summary-file
misc.nssummary=['_results.txt'];

%Run the nested sampling algorithm for all models and compile results
[results] = ns_processdataset(obs,models,misc);

%Save the results and the original data
%path=[misc.data_id,'_output'];
%save(path,'results')
%save(path,'data','-append')

