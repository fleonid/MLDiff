%%
%   smoothing zero-order spline for z, following NumDiff_v26 
%   created 2025-07-09
%
function [Q,C]=zero_order_spline_step_QC(h_K)
%
 K=length(h_K)+1;
%
 Q=diag([0;h_K]);
%
 C=zeros(K,K); 
 C(:,1)=ones(K,1);
 for j=2:K
     C(j:K,j)=h_K(j-1)*ones(K-j+1,1);
 end
% 
 return