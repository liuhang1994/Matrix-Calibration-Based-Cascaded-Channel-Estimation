function [NMSE]=NMSE_CAL(True,Estimate)
% Compute NMSE by eq. (40a)
NMSE=(norm(True-Estimate,'fro')/norm(True,'fro'))^2;

end