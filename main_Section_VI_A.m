clear;
clc;
warning('off');
basePath = [fileparts(mfilename('fullpath')) filesep];
addpath([basePath 'Replica_Library']);
addpath([basePath 'MP_Library']);
addpath([basePath 'Model_Generation_Library']);
libopt=[];
libopt.tau_N_inverse=Inf; %SNR level
libopt.M=51;
libopt.M_prime=64;                        %# of sampling grid for BS AoAs M'
libopt.T=60;                                   % pliot length
libopt.K=40;                                     %# of users
libopt.L1=4;                                     %vertical dimension of the LIS
libopt.L2=5;                                     %horizontal dimension of the LIS
libopt.L1_prime=4;
libopt.L2_prime=5;
libopt.L=libopt.L1*libopt.L2;          % # of antennas at the LIS
libopt.L_prime=libopt.L1_prime*libopt.L2_prime;
libopt.lambdaG=0.1;                    % lambda_u
libopt.lambdaS=0.05;
libopt.tauS=1;
libopt.tauG=1;
libopt.optiter=500; % I_max for message passing
libopt.trails=5000; % number of monte carlo trials
SNRlist=1;
libopt.tau_N_inverse=input('input the value of 1/tau_N\n');

libopt.pathstr=[basePath 'DATA/SNR_' num2str(libopt.tau_N_inverse) '.mat'];
Final_G=zeros(length(SNRlist),libopt.trails);
Final_S=zeros(length(SNRlist),libopt.trails);

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
%intial
system=[];
%noise power
system.nuw=10^(-libopt.tau_N_inverse/10);
for t=1:length(SNRlist)
    % libopt.SNR=SNRlist(t); %SNR level
    fprintf('tau_N_inverse: %d\n',libopt.tau_N_inverse)
    
    for i=1:libopt.trails
        if mod(i,1)==0
            fprintf('Iteration: %d\n',i)
        end
        %generate channel matrix
        %generate S
        G=S_GEN(L_prime,K,libopt.lambdaG)*sqrt(libopt.tauG);
        S=S_GEN(M_prime,L_prime,libopt.lambdaS)*sqrt(libopt.tauS);
        X=(randn(K,T)+1j*randn(K,T))/sqrt(2);
        H_RB_bar=(randn(M,L)+1j*randn(M,L))/sqrt(2);
        noise=sqrt(system.nuw/2)*(randn(M,T)+1i*randn(M,T));
        A_B=A_GEN(M,M_prime);
        varphi=-1+1/L1_prime:2/L1_prime:1-1/L1_prime;
        varsigma=-1+1/L2_prime:2/L2_prime:1-1/L2_prime;
        A_R=F_GEN(L1,L2,varphi,varsigma);
        
        
        system.S=S;
        system.G=G;
        system.X=X;
        system.A_B=A_B;
        system.R=A_R'*A_R;
        system.H_0=H_RB_bar*A_R;
        %     system.H_0=(randn(M,L_prime)+1j*randn(M,L_prime))/sqrt(2);
        Z=(system.H_0+system.A_B*system.S*system.R)*system.G;
        Q=Z*X;
        Y=Q+noise;
        system.Y=Y;
        
        
        
        Fin_OPT=MessagePassing(system,libopt);
        Final_G(t,i)=MSE_Compute(G,Fin_OPT.ghat);
        Final_S(t,i)=MSE_Compute(S,Fin_OPT.shat);
        if mod(i,1)==0
            fprintf('Updated:Su NMSE: %f, Sr NMSE: %f\n',10*log10(Final_G(t,i)),10*log10(Final_S(t,i)));
        end
        
        
    end
end
MSE_G_simulation=10*log10(mean(Final_G,2));
MSE_S_simulation=10*log10(mean(Final_S,2));
save(libopt.pathstr,'libopt','Final_G','Final_S','MSE_G_simulation','MSE_S_simulation')