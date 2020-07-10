function [estFin]= MessagePassing(system,libopt)
% Implementation of Algorithm 1
% Input: system contains the known realizations
%        libopt contains the system parameter
% Output:
% A structure that contains the estimates



% intialization

L_prime = libopt.L_prime;
K = libopt.K;
M_prime = libopt.M_prime;
M = libopt.M;
Y=system.Y;
A=system.A_B;
X=system.X;
H_0=system.H_0;
R=system.R;

% initalize a structure that contains the algorithm control parameters
opt = MPOpt;

% I_max
opt.nit = libopt.optiter;
% early stopping tolerance, epsilon
opt.tol=1e-4;


% see MP_Library/MPOpt.m
opt.adaptStep=1;
opt.pvarMin=1e-5;
opt.SvarMin=1e-5;
opt.nitMin=1;
% Initial Setup

% structure computing the prior related functions
gGbase = CAwgnEstimIn(0, libopt.tauG); % Complex Guassian Prior 
gSbase = CAwgnEstimIn(0, libopt.tauS); % Complex Guassian Prior 
gG = SparseScaEstim(gGbase,libopt.lambdaG); % Bernouli Complex Guassian Prior 
gS = SparseScaEstim(gSbase,libopt.lambdaS); % Bernouli Complex Guassian Prior 
gOut = CAwgnEstimOut(Y, system.nuw); % AWGN output prior 
 


% inital value for iteration 0
opt.ghat0=zeros(L_prime,K);
opt.shat0=S_GEN(M_prime,L_prime,libopt.lambdaS)*sqrt(libopt.tauS);
opt.what0=H_0+A*opt.shat0*R;
opt.zhat0=opt.what0*opt.ghat0;
opt.gvar0=ones(L_prime,K);
opt.svar0=ones(M_prime,L_prime);
opt.wvar0=ones(M,L_prime);
opt.zvar0=ones(M,K);

%message passing iteration
[estFin,~]=MessagePassing_iteration(R,X,H_0,A,gOut, gG, gS, opt);
