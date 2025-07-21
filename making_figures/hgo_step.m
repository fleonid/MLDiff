function X2=hgo_step(h_K,y_K,L1,ep)
%
% Levant's basic algorithm
%
K=length(y_K);
%
 L=L1+0.01;
 if ep==[],
   ep=1e-2;
 end
 k1=2/ep;
 k2=1/(ep^2);
%
 E=zeros(K,1);
 X1=zeros(K,1);
 X2=zeros(K,1);
%
nj=1;
for i=2:K
   E(i-1)=X1(i-1)-y_K(i-1);
   X1j1=X1(i-1);
   X2j1=X2(i-1);
   for j=1:nj 
     h=h_K(i-1)/nj;   
     X1j=X1j1 + h*( -k1*E(i-1) + X2j1 ); % try constant h instead of H_K
     X2j=X2j1 + h*( -k2*E(i-1) );
     X1j1=X1j;
     X2j1=X2j;
   end
   X1(i)=X1j;
   X2(i)=X2j;
end
%
for i=1:K
   X2(i)=min(max(X2(i),-L),L); % saturation
end
%
return