function z_K=z_from_p(p_K,h_K)
 K=length(p_K)+1;
   z_K=zeros(K,1); % i=1:K
%  z_K(1)=-p_K(1)/2; % z'=0
 for i=2:K-1
     z_K(i)=-1/2*(h_K(i)*p_K(i-1)+p_K(i)*h_K(i-1))/(h_K(i)+h_K(i-1));
 end
  z_K(1)=-z_K(2)-p_K(1); %  z''=0
  ddzK=0;
  z_K(K)=-z_K(K-1)-p_K(K-1)+h_K(K-1)^2*ddzK/2; %  z''=ddzK
return
