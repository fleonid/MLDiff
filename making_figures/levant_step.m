function X2=levant_step(h_K,y_K,L2)
%
% Levant's basic algorithm
%
K=length(y_K);
%
 L=L2+0.01;  % L=2*pi+1.5;
 k1=1.5*sqrt(L);
 k2=1.1*L;
%
 E=zeros(K,1);
 X1=zeros(K,1);
 X2=zeros(K,1);
%
 nj=1;
%
for i=2:K
   E(i-1)=X1(i-1)-y_K(i-1);
   sn=sign(E(i-1));
   X1j1=X1(i-1);
   X2j1=X2(i-1);
   for j=1:nj 
     h=h_K(i-1)/nj;   
     X1j=X1j1 + h*( -k1*sqrt(abs(E(i-1)))*sn + X2j1); 
     X2j=X2j1 + h*( -k2*sn );
     X1j1=X1j;
     X2j1=X2j;
   end
   X1(i)=X1j;
   X2(i)=X2j;
end
%
return