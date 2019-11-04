sigma=2;
mu=[0 0];
sigma_mn=0; %standard deviation for the measurement noise
N_steps=200;
dim=size(mu,2);
obs=sigma*randn(N_steps,dim)+mu;
data=[zeros(1,dim); cumsum(obs)];
data=data + sigma_mn*randn(N_steps+1,dim);  % Add measurement noise
save('trajectory_data.txt','data','-ascii');

