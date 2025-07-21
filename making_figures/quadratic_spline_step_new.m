function [x_1,p_K,z_K]=quadratic_spline_step_new(h_K,y_K,lambda)
%
K=length(y_K);
%
Y=y_K;
%
% find p_K
%
ZPmap=zeros(K,K); % [z1,z2,...,zK]=ZPmap*[p1,...,p_{K-1},zK]
%
for i=2:K-1
   ZPmap(i,i)=-1/2*(h_K(i-1))/(h_K(i)+h_K(i-1));
   ZPmap(i,i-1)=-1/2*(h_K(i))/(h_K(i)+h_K(i-1));
end
ZPmap(1,1)=-(h_K(1)+(h_K(2)/2))/(h_K(1)+h_K(2)); % z''=0
ZPmap(1,2)=(h_K(1)/2)/(h_K(1)+h_K(2)); % z''=0
%
%ZPmap(K,K-1)=-(1/2*h_K(K-2)+h_K(K-1))/(h_K(K-1)+h_K(K-2)); % z''=0
%ZPmap(K,K-2)=1/2*(h_K(K-1))/(h_K(K-1)+h_K(K-2)); % z''=0
ZPmap(K,K-1)=1;
%
 IPmap=zeros(K,K);
for j=1:K-1
   IPmap(j+1:K,j)=-h_K(j)/6*ones(K-j,1);
end
%
IZmap=zeros(K,K);
%
 for j=2:K
    IZmap(j:K,j-1)=IZmap(j:K,j-1)+h_K(j-1)/3; 
    IZmap(j:K,j)=IZmap(j:K,j)+h_K(j-1)/3;
 end
C=zeros(K,K+1);
C(:,1)=ones(K,1); 
C(:,2:K+1)=IPmap+IZmap*ZPmap; 
%
Z1Pmap=zeros(K-1,K-1);
for i=2:(K-1)
    Z1Pmap(i,i-1)=1/(h_K(i-1)+h_K(i));
    Z1Pmap(i,i)=-Z1Pmap(i,i-1);
end
 Z1Pmap(1,:)=Z1Pmap(2,:);
 B1=[zeros(K-1,1) Z1Pmap];
 Q0=B1'*diag(h_K)*B1;
%
 Z2Pmap=diag(2./h_K.^2)*(ZPmap(1:K-1,1:K-1)+ZPmap(2:K,1:K-1)+eye(K-1));
 B2=[zeros(K-1,1) Z2Pmap];
 Q2=B2'*diag(h_K.^3/3)*B2;
%
 Q1=(B2'*diag(h_K.^2)*B1+B1'*diag(h_K.^2)*B2)/2; % symmetric part
%
 Q=Q0+Q2+Q1;
%
% PQ=(Y'*C/(C'*C+lambda*Q))';
 PQ=(C'*C+lambda*Q)\C'*Y;
%
 x_1=PQ(1);
 p_K=PQ(2:K);
 z_K=PQ(K+1);
%
return