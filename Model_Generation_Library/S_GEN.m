function S=S_GEN(L,K,lambda)
% generate a L*K Bernouli complex Gaussian matrix with
% the variance of non-zero entris=1 and Bernouli parameter
% (sparsity)=lambda

%generate S follow iid standard normal
S=sqrt(1/2).*(randn(L,K)+1j*randn(L,K));
%dropout by probability (1-lambda)
if numel(lambda)==1
    S=S.*(rand(L,K) <= lambda);
elseif numel(lambda)==K
    S=S.*(rand(L,K) <= repmat(lambda,L,1));
end


end