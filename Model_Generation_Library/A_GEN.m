function [A,theta]=A_GEN(N,L)
% This function generates a N\times L BS sampling basis with a unfirom sampling grid

%uniform sampling
reslu=2/L;
theta=(-1+reslu/2:reslu:1-reslu/2);
%sterring vector
A = exp(-1i*pi*(0:N-1)'*theta)/sqrt(N);
end