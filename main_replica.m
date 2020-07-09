% This script can produce the solid pruple curve in Fig. 5
% Specifically, it provides an iterative algorithm to compute the
% asymptotic MSEs by computing the fixed-point of eq. (37)
% the result MSE_G_ana and MSE_S_ana are stored in DATA/VIA_Analytical.mat


% Note that different to eq. (37), here we allow L_prime \neq L
% All discussed in Section V-B, this leads to a loose bound
tau_N_inverse=-40:2:-20;
K=40;
M=round(1.28*K);
M_prime=round(1.6*K);
T=round(1.5*K);
L=round(K*0.5);
L_prime=round(K*0.5);
lambdaG=0.1;   
lambdaS=0.05;  
tauS=1;
tauG=1;
tauX=1;
tau_H_0=1;

[MSE_S_ana,MSE_G_ana]=replica_iteration(tau_N_inverse,...
    K,M,M_prime,T,L,L_prime,lambdaS,lambdaG,tauS,tauG,tauX,tau_H_0);

basePath = [fileparts(mfilename('fullpath')) filesep];
save([basePath 'DATA/VIA_Analytical.mat'],'tau_N_inverse','MSE_S_ana','MSE_G_ana')