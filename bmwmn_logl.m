function logl = bmwmn_logl(obs,mu_step,var_step,var_mn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim=size(mu_step,2);
mu=mu_step;
var=var_step+2*var_mn;
logl = -log(sqrt(2*pi*var))*dim-sum((obs(1,:)-mu).^2)/(2*var);

for t=2:length(obs)
   mu = mu_step-var_mn/var*(obs(t-1,:)-mu);
   var = var_step+var_mn*(2-var_mn/var);
   logl = logl-log(sqrt(2*pi*var))*dim-sum((obs(t,:)-mu).^2)/(2*var);
end

