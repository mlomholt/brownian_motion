%function bm_skel()  %Use this line if you do not want the script to interfere with your workspace
clear models misc % This avoids interference with previously defined versions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This skeleton script carries out Bayesian inference on a trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

%%% Load the trajectory data and convert from positions to steps
data=load('trajectory_data.txt');
obs=diff(data);

%%% There are two models:
% Model 1 is Brownian motion (theta(1)=sigma is the standard deviation for one step)
% Model 2 is Brownian motion with a drift (theta(2:3)=[mu_x mu_y] is the drift per step)

%%% Specify the log-likelihood functions for the models
log_normal=@(obs,mu,var) -log(sqrt(2*pi*var))*numel(obs)-sum(sum((obs-mu).^2))/(2*var);
models{1}.logl=@(obs,theta) log_normal(obs,0,theta(1)^2);
models{2}.logl=@(obs,theta) log_normal(obs,theta(2:3),theta(1)^2);

%%% Specify priors for the models, i.e., the inverse of the cumulative prior distributions
%
%% Prior for Model 1
inv_normal=@(u) sqrt(2)*erfinv(2*u-1);
span_orders_sigma=1; %width in orders of magnitude below or above typical_sigma
typical_sigma=1;
inv_sigma=@(u) typical_sigma*10^(span_orders_sigma*inv_normal(u));
models{1}.invprior=@(u) inv_sigma(u);
%
%% Prior for Model 2 (using Model 1's prior for sigma)
width_mu=10;
inv_mu=@(u) width_mu*inv_normal(u);
models{2}.invprior=@(u) [inv_sigma(u(1)), inv_mu(u(2)), inv_mu(u(3))];

%%% Specify the random u-generators (that provide inputs for the invprior-functions)
models{1}.genu=@() rand(1,1);
models{2}.genu=@() rand(1,3);

%%% Labels for the parameters in the summary text-file
models{1}.labels={'sigma:'};
models{2}.labels={'sigma:','mu_x: ','mu_y: '};

%%% The filename for a text summary of the nested sampling analysis
misc.nssummary=['bm_results.txt'];

%%% Run the nested sampling algorithm for all models and compile results
[results] = ns_main(obs,models,misc);

