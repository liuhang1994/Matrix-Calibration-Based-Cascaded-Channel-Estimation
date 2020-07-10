% This script produces the data for the purple dashed curves in Fig. 5 for
% given noise power (tau_N);
% The results are MSE_G_simulation and MSE_S_simulation, which is also
% stored in DATA/VIA_Simulation.mat;
% To save the running time, here we can set the number of Monte Carlo
% trials as a small number.
% To fully recover the plots in Fig. 5, one can change libopt.trails
% to 5000;
clear;
clc;
warning('off');
%% add script path
basePath = [fileparts(mfilename('fullpath')) filesep];
addpath([basePath 'Replica_Library']);
addpath([basePath 'MP_Library']);
addpath([basePath 'Model_Generation_Library']);

%% libopt contains system parameter that will be passed into message passing functions
libopt=[];

% change to 5000 for smooth result
libopt.trails=5000;    % number of monte carlo trials



libopt.M=51;        % number of antennas at BS
libopt.M_prime=64;  % dimension of sampling grid at BS
libopt.T=60;        % pliot length
libopt.K=40;        % # of users
libopt.L1=4;        % vertical dimension of the LIS
libopt.L2=5;        % horizontal dimension of the LIS
libopt.L1_prime=4;  % vertical dimension of sampling grid at RIS
libopt.L2_prime=5;  % horizontal dimension of sampling grid at RIS

libopt.L=libopt.L1*libopt.L2;          % # of antennas at the LIS
libopt.L_prime=libopt.L1_prime*libopt.L2_prime;  % total sampling dimension
libopt.lambdaG=0.1;                    % lambda_G defined in eq.(13), sparsity of G
libopt.lambdaS=0.05;                   % lambda_S defined in eq.(12), sparsity of S
libopt.tauS=1;                         % tau_S, nonzero variance
libopt.tauG=1;                         % tau_G, nonzero variance
libopt.optiter=500; % I_max, maximum allowable number of algorithm 1


% input noise power by user
% For example: the range in Fig. 5 is -40:2:-20
SNRlist=input('input the range of 1/tau_N (in the form of list)\n');




% rename the parameters for convenience
K=libopt.K;
T=libopt.T;
M=libopt.M;
M_prime=libopt.M_prime;
L_prime=libopt.L_prime;
L=libopt.L;
L1_prime=libopt.L1_prime;
L1=libopt.L1;
L2_prime=libopt.L2_prime;
L2=libopt.L2;


%system contains the realization of each Monte Carlo trail
system=[];
% store MSEs of G and S for each trial
Final_G=zeros(length(SNRlist),libopt.trails);
Final_S=zeros(length(SNRlist),libopt.trails);
libopt.pathstr=[basePath 'DATA/VIA_Simulation.mat'];
for t=1:length(SNRlist)
    libopt.tau_N_inverse=SNRlist(t);
    %file name
    
    fprintf('tau_N_inverse: %d\n',libopt.tau_N_inverse)
    %noise power
    system.nuw=10^(-libopt.tau_N_inverse/10);
    for i=1:libopt.trails
        
        %verbose print
        if mod(i,1)==0
            fprintf('Iteration: %d\n',i)
        end
        
        %generate channel matrices
        G=S_GEN(L_prime,K,libopt.lambdaG)*sqrt(libopt.tauG); %eq. (13)
        S=S_GEN(M_prime,L_prime,libopt.lambdaS)*sqrt(libopt.tauS); %eq. (12)
        % Angular basis with grids uniformly sample [-1,1]
        A_B=A_GEN(M,M_prime); %A_B, (over-complete) DFT
        %A_R, (over-complete) 2D-DFT
        A_R=F_GEN(L1,L2,...
            -1+1/L1_prime:2/L1_prime:1-1/L1_prime,-1+1/L2_prime:2/L2_prime:1-1/L2_prime);
        X=(randn(K,T)+1j*randn(K,T))/sqrt(2);
        
        system.S=S; % S in eq. (10)
        system.G=G; % G in eq. (10)
        system.X=X; % X in eq. (10)
        system.A_B=A_B; % A_B in eq. (10)
        system.R=A_R'*A_R; % R in eq. (10)
        
        % H_RB_bar=(randn(M,L)+1j*randn(M,L))/sqrt(2);
        % system.H_0=H_RB_bar*A_R;
        system.H_0=(randn(M,L_prime)+1j*randn(M,L_prime))/sqrt(2); %H_0 in eq.(10)
        % Note: distribution of H_0 has no significant influence on the
        % final result, as long as its entries are i.i.d.
        
        
        %received signal
        noise=sqrt(system.nuw/2)*(randn(M,T)+1i*randn(M,T));
        Z=(system.H_0+system.A_B*system.S*system.R)*system.G;
        Q=Z*X;
        Y=Q+noise;
        system.Y=Y;
        
        
        % Run Algorithm 1
        Fin_OPT=MessagePassing(system,libopt);
        % Compute MSEs
        Final_G(t,i)=MSE_Compute(G,Fin_OPT.ghat);
        Final_S(t,i)=MSE_Compute(S,Fin_OPT.shat);
        
        % verbose print
        if mod(i,1)==0
            fprintf('Updated:Su NMSE: %f, Sr NMSE: %f\n',...
                10*log10(Final_G(t,i)),10*log10(Final_S(t,i)));
        end
        
        
    end
end
% Average MSEs
MSE_G_simulation=10*log10(mean(Final_G,2));
MSE_S_simulation=10*log10(mean(Final_S,2));
% save the data
save(libopt.pathstr,'libopt','Final_G','Final_S','MSE_G_simulation','MSE_S_simulation')

