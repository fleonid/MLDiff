% FILE: make_data.m
% DESCRIPTION:
%   Generates synthetic signal data with non-uniform sampling and Gaussian noise.
%
% INPUTS:
%   h         - average sampling step
%   dh        - relative variation in sampling step
%   t_initial - initial time
%   t_final   - final time
%   sigma   - standard deviation of Gaussian noise
%
% OUTPUTS:
%   K      - number of samples
%   t_K    - sample time points
%   h_K    - actual sampling intervals
%   x_K    - true signal values at t_K
%   dx_K   - true derivative values at t_K
%   y_K    - noisy measurements
%   L0num  - max absolute value of signal
%   L1num  - max absolute value of derivative
%   L2num  - max absolute value of second derivative
%   t      - fine time grid
%   x      - true signal on fine grid
%   dx     - true derivative on fine grid
%
% AUTHOR:
%   Leonid B. Freidovich
%
% LAST UPDATED:
%   2025-07-18
function [K,t_K,h_K,x_K,dx_K,y_K,L0num,L1num,L2num,t,x,dx]=make_data(h,dh,t_initial,t_final,sigma)
%
h_max=(1+dh)*h;
h_min=(1-dh)*h;
t_K(1)=t_initial;
i=0;
while 1
i=i+1;
    h_K(i)=h_min+(h_max-h_min)*rand(); % uniform noise
    t_K(i+1)=t_K(i)+h_K(i);
    if t_K(i+1)>t_final
        K=i+1;
   %     t_K(K)=t_final;
   %     h_K(K)=t_K(K)-t_K(K-1);
        break
    end
end
%
 t_K=t_K'; % making column from row
 h_K=h_K';
%
%%
% figure(1), subplot(211), plot(T_K), grid, title('time samples'), subplot(212), histogram(T_K,20);
%%
x_K=zeros(K,1);
y_K=zeros(K,1);
dx_K=zeros(K,1);
n_K=zeros(K,1);
%
for i=1:K
  [x_K(i),dx_K(i)]=x_fun(t_K(i));
  n_K(i)=sigma*randn(); 
  y_K(i)=x_K(i)+n_K(i);
end
%
t=(t_K(1):h_min/10:t_K(K))';
x=zeros(size(t)); dx=zeros(size(t));
L2num=0;
for i=1:length(t)
  [x(i),dx(i)]=x_fun(t(i));
  if i>1 
      L2numprev=L2num;
      temp=abs(dx(i)-dx(i-1))/abs(t(i)-t(i-1));
      L2num=max(temp,L2numprev);
  end
end
L0num=max(abs(x));
L1num=max(abs(dx));
return
