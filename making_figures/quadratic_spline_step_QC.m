function [Q,C]=quadratic_spline_step_QC(h_K)
%
 K=length(h_K)+1;
%
  if K<4
     disp(' '), fprintf('K = %i', K),  disp(' '), 
     disp(' ERROR:  K must be bigger than 3 '), 
     disp(' '),
 end
%
% find p_K
%
ZPmap=zeros(K,K-1);
%
for i=2:K-1
   ZPmap(i,i)=-1/2*(h_K(i-1))/(h_K(i)+h_K(i-1));
   ZPmap(i,i-1)=-1/2*(h_K(i))/(h_K(i)+h_K(i-1));
end
%ZPmap(K,K-1)=-1/2;  % z'=0
%ZPmap(1,1)=-1/2; % z'=0
ZPmap(1,1)=-(h_K(1)+(h_K(2)/2))/(h_K(1)+h_K(2)); % z''=0
ZPmap(1,2)=(h_K(1)/2)/(h_K(1)+h_K(2)); % z''=0
%
ZPmap(K,K-1)=-(1/2*h_K(K-2)+h_K(K-1))/(h_K(K-1)+h_K(K-2)); % z''=0
ZPmap(K,K-2)=1/2*(h_K(K-1))/(h_K(K-1)+h_K(K-2));
%
 IPmap=zeros(K,K-1);
for j=1:K-1
   IPmap(j+1:K,j)=-h_K(j)/6*ones(K-j,1);
end
%
IZmap=zeros(K,K);
%for j=1:K-1
%   IZmap(j+1:K,j)=h_K(j)/3*ones(K-j,1);
%end
%
%IZmap(1:K,2:K)=IZmap(1:K,2:K)+IZmap(1:K,1:K-1);
%
 for j=2:K
    IZmap(j:K,j-1)=IZmap(j:K,j-1)+h_K(j-1)/3;  % corrected on 2025-03-18
    IZmap(j:K,j)=IZmap(j:K,j)+h_K(j-1)/3;
 end
C=zeros(K,K);
C(:,1)=ones(K,1); 
C(:,2:K)=IPmap+IZmap*ZPmap; 
%
% B=zeros(K-1,K);
% for i=1:(K-2)
%     B(i,i+1)=1/(h_K(i)+h_K(i+1));
%     B(i,i+2)=-B(i,i+1);
% end
%     B(K-1,K-1)=1/(h_K(K-2)+h_K(K-1)); % z''=0
%     B(K-1,K)=-B(K-1,K-1);
%
Z1Pmap=zeros(K-1,K-1);
for i=2:(K-1)
    Z1Pmap(i,i-1)=1/(h_K(i-1)+h_K(i));
    Z1Pmap(i,i)=-Z1Pmap(i,i-1);
end
   Z1Pmap(1,:)=Z1Pmap(2,:);
%    B(K-1,K-1)=sym(1)/(h(K-2)+h(K-1)); % z''=0
%    B(K-1,K)=-B(K-1,K-1);
 B=[zeros(K-1,1) Z1Pmap];
 Q0=B'*diag(h_K)*B;
%
 Z2Pmap=diag(2./h_K.^2)*(ZPmap(1:K-1,1:K-1)+ZPmap(2:K,1:K-1)+eye(K-1));
 B2=[zeros(K-1,1) Z2Pmap];
 Q2=B2'*diag(h_K.^3/3)*B2;
%
 Q1=(B2'*diag(h_K.^2)*B+B'*diag(h_K.^2)*B2)/2; % symmetric part
%
 Q=Q0+Q2+Q1;
%
% PQ=(Y'*C/(C'*C+lambda*Q))';
% PQ=(C'*C+lambda*Q)\C'*Y;
%
%x_1=PQ(1);
%p_K=PQ(2:K);
%
return