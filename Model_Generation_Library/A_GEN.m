function [A,theta_f]=A_GEN(N,L)
%     L=libopt.L;
%     N=libopt.N;
    %uniform sampling
    reslu=2/L;
    theta=(-1+reslu/2:reslu:1-reslu/2);
    %add noise
%     beta=(2*rand(1,L)-1)/L;
    theta_f=theta;
    %sterring vector
    A = exp(-1i*pi*(0:N-1)'*theta_f)/sqrt(N);
end