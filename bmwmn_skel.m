%function bmwmn_skel()  %Use this line if you do not want the script to interfere with your workspace
clear models misc % This avoids interference with previously defined models arrays

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This skeleton script carries out Bayesian inference on a trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

%%% Load the trajectory data and convert from positions to steps
data=load('trajectory_data.txt');
obs=diff(data);

%%% There are three models:
% Model 1 is Brownian motion (theta(1)=sigma is the standard deviation for one step)
% Model 2 is Brownian motion with a drift (theta(2:3)=[mu_x mu_y] is the drift per step)
% Model 3 is Brownian motion with measurement error (theta(2) is the standard deviation for the measurement noise).
% Model 4 is Brownian motion with drift and measurement error (theta(2:3)=[mu_x mu_y]).

%%% Specify the log-likelihood functions for model 1 and 2
log_normal=@(obs,mu,var) -log(sqrt(2*pi*var))*numel(obs)-sum(sum((obs-mu).^2))/(2*var);
models{1}.logl=@(obs,theta) log_normal(obs,[0,0],theta(1)^2);
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

%% Prior and log-likelihood for model 3
span_orders_mn=1; %width in orders of magnitude below or above typical_mn
typical_mn=1;
inv_mn=@(u) typical_mn*10^(span_orders_mn*inv_normal(u));
models{3}.invprior=@(u) [inv_sigma(u(1)), inv_mn(u(2))];
models{3}.logl=@(obs,theta) bmwmn_logl(obs,[0, 0],theta(1)^2,theta(2)^2);

%% Everything for for model 4
models{4}.invprior=@(u) [inv_sigma(u(1)), inv_mu(u(2)), inv_mu(u(3)), inv_mn(u(4))];
models{4}.logl=@(obs,theta) bmwmn_logl(obs,theta(2:3),theta(1)^2,theta(4)^2);

%%% Specify the random u-generators (that provide inputs for the invprior-functions)
models{1}.genu=@() rand(1,1);
models{2}.genu=@() rand(1,3);
models{3}.genu=@() rand(1,2);
models{4}.genu=@() rand(1,4);

%%% Labels for the parameters in the summary text-file
models{1}.labels={'sigma:'};
models{2}.labels={'sigma:','mu_x: ','mu_y: '};
models{3}.labels={'sigma:','mn:'};
models{4}.labels={'sigma:','mu_x: ','mu_y: ','mn:'};

%%% The filename for a text summary of the nested sampling analysis
misc.nssummary=['bmwmn_results.txt'];

%%% Run the nested sampling algorithm for all models and compile results
[results] = ns_main(obs,models,misc);

