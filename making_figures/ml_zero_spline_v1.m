% FILE: ml_zero_spline_v1.m
% DESCRIPTION:
%   Implements the zero-order spline-based differentiator and compares it with Levant's method.
%
% CONTEXT:
%   Based on the method described in:
%   "Numerical Differentiation under Coarse Non-uniform Sampling and Gaussian Noise"
%   by K. E. Avrachenkov and L. B. Freidovich, submitted to IEEE TAC, 2025.
%
% INPUTS:
%   None directly; parameters such as epsilon, h, t_initial, t_final are set in the script.
%
% OUTPUTS:
%   Plots comparing analytical and estimated derivatives using zero-order spline and Levant's method.
%
% AUTHOR:
%   Leonid B. Freidovich (leonid.freidovich@umu.se)
%
% LAST UPDATED:
%   2025-07-18
 clearvars
 close all
%
% diary ml_quadratic_v13.txt
%
 disp(' '), disp(datetime('now')),
%
%  epsilon=1e-5; % amplitude of noise0
%  epsilon=0.000001; % amplitude of noise0
epsilon=0.001; 
%disp('amplitude of noise'), disp(epsilon)
disp(' '), fprintf('amplitude of noise  epsilon = %f', epsilon), disp(' '),
%
%h=0.05; % average sampling step
%h=0.01; % average sampling step
h=0.05;
%  h=0.001;
disp(' '), fprintf('average sampling  h = %f', h), disp(' '),
%disp('average sampling'), disp(h)
%
 t_initial=0;
%  t_final=0.05;
%
% t_final=2.13; % final time
%   t_final=1.22; % final time
% t_final=0.3; % final time
% t_final=25;
% t_final=10.2;
% t_final=9.75;
  t_final=5.11; 
% t_final=2.5;
%  t_final=150;

%disp('final time'), disp(t_final)
%disp(' '), fprintf('final time,  t_final = %f', t_final), disp(' '),
%
dh=0.5; % 50% tikme step variation
%
% make data
%
 [K,t_K,h_K,x_K,dx_K,y_K,L0num,L1num,L2num,t,x,dx]=make_data(h,dh,t_initial,t_final,epsilon);
%
disp(' '), fprintf('t_K(1)=%f, x_K(1)=%f, dx_K(1)=%f',t_K(1), x_K(1), dx_K(1)), disp(' '),
disp(' '), fprintf('t_K(K)=%f, x_K(K)=%f, dx_K(K)=%f',t_K(K), x_K(K), dx_K(K)), disp(' '),
% figure(2), 
% subplot(211), plot(t_K,x_K,'x',t_K,y_K,'o',t,x), 
% grid, legend('x_k=x(t_k)', 'y_k=x(t_k)+n(t_k)', 'x=x(t)','Location','northwest')
% subplot(212), plot(t_K,dx_K,'x',t_K,n_K,'*',t,dx), 
% grid, legend('dx_k=dx/dt(t_k)', 'n_k=n(t_k)', 'dx=dx(t)','Location','southeast')
%
 lambda=1e-4; % empirically found
% lambda=2.2e-7; % theoretically optimal
% lambda=1e-6;
%
 disp(' '), fprintf('lambda = %g', lambda), disp(' ')
%
 C=zeros(K,K); Q=diag([0;h_K]);
 C(:,1)=ones(K,1);
 for j=2:K
     C(j:K,j)=h_K(j-1)*ones(K-j+1,1);
 end
% 
 A=C'*C+lambda*Q; 
 PQ=A\(C'*y_K); 
 x0=PQ(1);  
 z_K=PQ(2:K);
%
 X2=levant_step(h_K,y_K,L2num);
%
z=zeros(size(t));
%
 i=1; % start from interval 1: (t_1,t_2)
%
 for j=1:length(t)
     z(j)=z_K(i);
     if t(j)>t_K(i+1)
        i=i+1;
     end
 end
%
 figure,
 plot(t,dx,'r',t_K,dx_K,'r.',...
    t_K,X2,'k', t_K,X2,'k.',...
   t,z,'g'), grid,
 legend('$dx/dt$','$dx/dt$ at $t_k$',...
    'Levant diff.','Levant diff. at $t_k$',...
    'zero spline',...
    'Location','south','Interpreter','latex');
% disp(' '), disp('Variances:'), disp(' ')
% 
% t_tr=floor(length(X2)/4);
% disp(' - Levant'), disp(sqrt(mean((X2(t_tr:end)-dx_K(t_tr:end)).^2))),
% disp(' - Recursive Quadratic Spline'), disp(sqrt(mean((z_K(t_tr:end)-dx_K(t_tr:end)).^2))),
% disp(' - Quadratic Spline'), disp(sqrt(mean((z_Kp(t_tr:end)-dx_K(t_tr:end)).^2))),

return
%%
 clearvars;
 syms lambda
 assume(lambda, {'real','positive'})
%
 K=3;
 h = sym('h', [K-1 1]);
 assume(h, {'real','positive'})
%
 CK=sym(zeros(K,K)); QK=diag([sym(0),transpose(h)]);
 CK(:,1)=sym(ones(K,1));
 for j=2:K
     CK(j:K,j)=h(j-1)*sym(ones(K-j+1,1));
 end
% 
 AK=transpose(CK)*CK+lambda*QK; 
% 
 K=K+1;
 h = sym('h', [K-1 1]);
 assume(h, {'real','positive'})
%
 CKp=sym(zeros(K,K)); QKp=diag([sym(0),transpose(h)]);
 CKp(:,1)=sym(ones(K,1));
 for j=2:K
     CKp(j:K,j)=h(j-1)*sym(ones(K-j+1,1));
 end
% 
 AKp=transpose(CKp)*CKp+lambda*QKp; 
%
 Delta=AKp-[AK,sym(zeros(K-1,1));sym(zeros(1,K-1)),sym(0)]
%
 Delta_new=[1;h]*transpose([1;h])+[sym(zeros(K-1,K-1)),sym(zeros(K-1,1));sym(zeros(1,K-1)),lambda*h(K-1)]
 simplify(Delta-Delta_new)
 f=sym(1);
 if K==3
  f=h(1)+2*lambda;
 end
 if K==4 
  f=h(1)*h(2) + 2*(h(1)+h(2))*lambda+3*lambda^2;
 end
 if K==5
  f=h(1)*h(2)*h(3) + 2*(h(1)*h(2)+h(1)*h(3)+h(2)*h(3))*lambda +(3*h(1)+4*h(2)+3*h(3))*lambda^2+4*lambda^3;
 end
 collect(inv(AK)*f,lambda)
 %
 y = sym('y', [K 1]);
 assume(y, {'real'})
 syms x0
 assume(x0, {'real'})
 DK=y(1:K-1); DKp=y;
 %
 ZK=AK\(CK*DK);
 ZKp=AKp\(CKp*DKp);
 collect(ZKp(1)-ZK(1),y)
 collect(ZKp(2:end-1)-ZK(end),y)
 collect(ZKp(end),y)
