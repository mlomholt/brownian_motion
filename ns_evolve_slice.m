function [walker_new,new_step_mod]=ns_evolve_slice(obs,model,logLstar,walker,step_mod,mode)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses the slice sampling technique to find a new sample
% uniformly distributed inside
% the region of parameter space with higher likelihood than the minimum 
% requirement (logLstar). 
%
% Some of the arguments of the function are
% 
% obs - a 2xT matrix of observations
% walker - the walker that constitutes the starting point of the
%   MCMC process.
% logLstar - the minimum requirement for the likelihood of the new sample
% step_mod - a scalar value that regulates the average step length of the
%   subsequent call of the function. When the remaning parameter space
%   becomes small, the MCMC steps are adjusted in length to ensure a
%   success rate of the MCMC steps of about 50%.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize step_mode if run for the first time
if step_mod==0
  step_mod=ones(size(walker.u));
end
new_step_mod=step_mod;

% Attempt to generate new walker through independent guess
walker_new.u = model.genu();
walker_new.theta = model.invprior(walker_new.u);
walker_new.logl=model.logl(obs,walker_new.theta); % Calculates new likelihood 

% If it failed, run slice sampling
if ~isequal(mode,'ns')
  logLstar=log(rand)+walker.logl;
end
if(walker_new.logl <= logLstar)	% Do MCMC if likelihood requirement failed
   % Find a new walker via a MCMC process.
   % Return input values and logLstar if no steps succeed
   [walker_new.u,walker_new.theta,walker_new.logl,new_step_mod] = ns_slice(obs,walker.u,walker.logl,model.invprior,model.logl,logLstar,step_mod,new_step_mod,model.options.nsteps,mode,model.options.nsteps);
end

