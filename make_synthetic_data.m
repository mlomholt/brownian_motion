sigma=2;
mu=[0 0];
N_steps=200;
dim=size(mu,2);
obs=sigma*randn(N_steps,dim)+mu;
data=[zeros(1,dim); cumsum(obs)];
save('trajectory_data.txt','data','-ascii');

