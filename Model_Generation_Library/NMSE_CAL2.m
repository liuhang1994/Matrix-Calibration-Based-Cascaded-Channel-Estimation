function [NMSE]=NMSE_CAL2(True,Estimate)
% Compute NMSE by eq. (40b)

a1=abs(True-Estimate).^2;
a2=abs(True).^2;
NMSE=mean(sum(a1)./sum(a2));
end