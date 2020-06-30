function [rijN]=L1Problem2(rij,Constraints)
%L1  reweighted L1 minimization for Problem2
% This function uses reweighted l_1 minimization to minimize the total
% number of non-zero link allocations. 
%
% Inputs:
% - resource allocation (rij)
%
% Outputs:
% - reweighted resource allocation (rijN)

% Vera Somers, March 2020

eps=1E-6;
temp=1;

for i=1:50

rijO=double(rij);

optimize(Constraints,sum(sum(rij./(rijO+eps))),sdpsettings('verbose',0));
if abs(temp - sum(sum(double(rij)./(rijO+eps)))) <= 1E-5
    break
else
    
temp=sum(sum(double(rij)./(rijO+eps)));

end

rijN=double(rij);

end
