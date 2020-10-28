%function corr_skel()  %Use this line if you do not want the script to interfere with your workspace
clear models misc % This avoids interference with previously defined versions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This skeleton script carries out Bayesian inference on a trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

%%% Load the trajectory data and convert from positions to steps
data=load('trajectory_data.txt');
obs=diff(data);

%%% There are two models:
% Model 1 is Brownian motion (theta(1)=sigma is the standard deviation for one step)
% Model 2 is Brownian motion with different diffusion constants in different directions (theta(1:2)=[sigma_x sigma_y] is the standard deviations in each direction)
% Model 3 is Model 2 but including correlation (theta(1:3)=[sigma_x sigma_y rho] where rho is the correlation coefficient)

%%% Specify the log-likelihood functions for the models
log_normal=@(obs,mu,var) -log(2*pi*var)*numel(obs)/2-sum(sum((obs-mu).^2))/(2*var);
models{1}.logl=@(obs,theta) log_normal(obs,0,theta(1)^2);
log_normal_2D=@(obs,theta) -log(2*pi*theta(1)*theta(2)*sqrt(1-theta(3)^2))*size(obs,1)-sum(obs(:,1).^2/theta(1)^2+obs(:,2).^2/theta(2)^2-2*theta(3)*obs(:,1).*obs(:,2)/theta(1)/theta(2))/2/(1-theta(3)^2);
models{2}.logl=@(obs,theta) log_normal_2D(obs,[theta 0]);
models{3}.logl=@(obs,theta) log_normal_2D(obs,theta);

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
models{2}.invprior=@(u) [inv_sigma(u(1)), inv_sigma(u(2))];

%% Prior for Model 3 (using Model 1's prior for sigma)
models{3}.invprior=@(u) [inv_sigma(u(1)), inv_sigma(u(2)), 2*u(3)-1];

%%% Specify the random u-generators (that provide inputs for the invprior-functions)
models{1}.genu=@() rand(1,1);
models{2}.genu=@() rand(1,2);
models{3}.genu=@() rand(1,3);

%%% Labels for the parameters in the summary text-file
models{1}.labels={'sigma_x:'};
models{2}.labels={'sigma_x:','sigma_y:'};
models{3}.labels={'sigma_x:','sigma_y:','rho:'};

%%% The filename for a text summary of the nested sampling analysis
misc.nssummary=['corr_results.txt'];

%%% Run the nested sampling algorithm for all models and compile results
[results] = ns_main(obs,models,misc);

