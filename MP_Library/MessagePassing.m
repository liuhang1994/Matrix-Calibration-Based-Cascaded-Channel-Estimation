function [estFin]= MessagePassing(system,libopt)
%% intialization
L_prime = libopt.L_prime;
K = libopt.K;
M_prime = libopt.M_prime;
M = libopt.M;
Y=system.Y;
A=system.A_B;
X=system.X;
H_0=system.H_0;
R=system.R;

opt = MPOpt;
opt.nit = libopt.optiter;
opt.adaptStep=1;
opt.tol=1e-4;
opt.pvarMin=1e-5;
opt.SvarMin=1e-5;
opt.XvarMin=1e-5;
opt.nitMin=1;
%% Initial Setup
gGbase = CAwgnEstimIn(0, libopt.tauG);
gSbase = CAwgnEstimIn(0, libopt.tauS);
gG = SparseScaEstim(gGbase,libopt.lambdaG);
gS = SparseScaEstim(gSbase,libopt.lambdaS);
gOut = CAwgnEstimOut(Y, system.nuw);

opt.ghat0=zeros(L_prime,K);
opt.shat0=S_GEN(M_prime,L_prime,libopt.lambdaS)*sqrt(libopt.tauS);
opt.what0=H_0+A*opt.shat0*R;
opt.zhat0=opt.what0*opt.ghat0;

opt.gvar0=ones(L_prime,K);
opt.svar0=ones(M_prime,L_prime);
opt.wvar0=ones(M,L_prime);
opt.zvar0=ones(M,K);


%message passing
[estFin,~]=MessagePassing_iteration(R,X,H_0,A,gOut, gG, gS, opt);
%     n1=MSE_Compute(system.S_u,estFin3.uhat);
%     n2=MSE_Compute(system.S_r,estFin3.rhat);
%     fprintf('Su NMSE: %f, Sr NMSE: %f\n',10*log10(n1),10*log10(n2));

