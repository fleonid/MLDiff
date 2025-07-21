function [dx,dP,pK,Qp,Cp,invAp]=update_quadratic_spline(x0,P,h,yK,yKp,lambda,Q,C,invA)
%
%profile on
  K=length(h); % 
% [Q,C]=quadratic_spline_step_QC(h(1:K-1)); % computed recursively
 [dQ,dC]=updateQC(h(1:K));
%
% Compute Delta
% Delta = dC'*dC+dC'*[[C(K,1:K),0];zeros(1,K+1)]+[[C(K,1:K)';0],zeros(K+1,1)]*dC+...
%        lambda*[zeros(K-2,3);eye(3,3)]*dQ*[zeros(3,K-2),eye(3,3)];
%
temp=[[C(K,1:K),0];zeros(1,K+1)]; 
Delta = dC'*(dC+temp)+temp'*dC; 
Delta(K-1:K+1,K-1:K+1)=Delta(K-1:K+1,K-1:K+1)+lambda*dQ;

% Compute Qp and Cp
%Qp = [Q, zeros(K, 1); zeros(1, K+1)];
Qp=zeros(K+1,K+1); Qp(1:K,1:K)=Q; Qp(K-1:K+1,K-1:K+1)=Qp(K-1:K+1,K-1:K+1)+dQ;

%Cp = [C, zeros(K, 1); zeros(1, K+1)];
Cp=zeros(K+1,K+1); Cp(1:K,1:K)=C; Cp(K:K+1,1:K+1)=Cp(K:K+1,1:K+1)+dC;

% Extract components of Delta=[...,Us;Vs,Cs] with Vs=Us'.
Us = Delta(1:K,K+1);
%Vs = Delta(K+1,1:K);
Cs = Delta(K+1,K+1);

% Efficient updates of invA using 3 times Sherman-Morrison-Woodbury
V1 = [dC(1,1:K)+C(K,1:K); dC(2,1:K)];
temp = invA*dC(:,1:K)';
invA1 = invA-temp/(eye(2,2)+V1*temp)*(V1*invA);
%
temp = invA1*[C(K,1:K)',zeros(K,1)];
invA2 = invA1-temp/(eye(2,2)+dC(1:2,1:K)*temp)*(dC(1:2,1:K)*invA1);
%
V3 = [zeros(2,K-2),eye(2,2)];
invAs=invA2-invA2(:,K-1:K)/(inv(dQ(1:2,1:2))/lambda+V3*invA2(:,K-1:K))*(V3*invA2);
%
invAsUs = invAs*Us; invAsUsinvAsVs=invAsUs*invAsUs';
invAn = invAs-(1/(-Cs+Us'*invAsUs))*(invAsUsinvAsVs);
%
% Compute update term
B = dC' * [yK; yKp] - Delta * [x0; P; 0];
upd = invAn * (B(1:K) - Us * B(K+1) / Cs);
%
% Extract updates
dx = upd(1);
dP = upd(2:K);
pK = (B(K+1) - Us' * upd) / Cs;

% Compute invAp 
invCc = 1/(Cs-Us'*invAsUs);
%invAp=zeros(K+1,K+1);
%invAp1 = [invAs + invCc*(((invAs * Us)*Vs)*invAs), -invCc*(invAs*Us); ...
%         -invCc*(Vs*invAs), invCc];
%
%invAp = [invAs+invCc*invAsUsinvAsVs, -invCc*invAsUs; -invCc*invAsUs', invCc];
invAp=zeros(K+1,K+1); 
invAp(1:K,1:K)=invAs+invCc*invAsUsinvAsVs; invAp(1:K,K+1)=-invCc*invAsUs; 
invAp(K+1,1:K)=invAp(1:K,K+1)'; invAp(K+1,K+1)=invCc;
%profile off
return