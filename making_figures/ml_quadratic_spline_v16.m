%%
%   Smoothing quadratic-spline-based differentiation, following NumDiff_v27
%   ["Numerical Differentiation under Coarse Non-uniform Sampling and Gaussian Noise"
%   by Konstantin E. Avrachenkov and Leonid B. Freidovich
%   (leonid.freidovich@umu.se)] 
%
%   updated 2025-07-18
%
 clearvars
 close all
%
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultAxesLineWidth', 1);  % Set default axes line width to 2 points
%
% diary ml_quadratic_v16.txt
 disp(' '), disp(datetime('now')),
%
% standard deviation of Gaussian noise
%
% sigma=0.0001; % Fig. 1a, 1b, 3
 sigma=1e-7; % Fig. 2a
% sigma=0.01; % Fig. 2b
% sigma=0.001; 
 disp(' '), fprintf('standard deviation of Gaussian noise sigma = %f', sigma), disp(' '),
%
%    average sampling step
%
% h=0.05; % Fig. 3
% h=0.01; % Fig. 1a
 h=0.001;  % Fig. 1b, 2a, 2b
% h=0.075;
disp(' '), fprintf('average sampling  h = %f', h), disp(' '),
%
 t_initial=0;
%  t_initial=0.05;
%
% t_final=2.13; % final time
%   t_final=1.22; % final time
% t_final=0.3; % final time
% t_final=25;
% t_final=10.2;
% t_final=9.75;
%  t_final=5.11; 
 t_final=1.95;
%  t_final=150;

disp(' '), fprintf('final time,  t_final = %f', t_final), disp(' '),
%
dh=0.5; % 50% time step variation
%
% make data
%
 [K,t_K,h_K,x_K,dx_K,y_K,L0num,L1num,L2num,t,x,dx]=make_data(h,dh,t_initial,t_final,sigma);
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
%  QUADRATIC SPLINES
%
% [x0,p_K]=quadratic_spline_step(h_K,y_K,lambda);
 [Q,C]=quadratic_spline_step_QC(h_K); 
 PQ=(C'*C+lambda*Q)\C'*y_K; 
 x0=PQ(1); 
 p_K=PQ(2:K);
 z_K=z_from_p(p_K,h_K);
%
 z=interp_z_from_z_K(t,t_K,h_K,z_K,p_K);
%
%
%  ZERO ORDER SPLINES
%
 [Q,C]=zero_order_spline_step_QC(h_K); 
 PQ=(C'*C+lambda*Q)\C'*y_K; 
 x0c=PQ(1); 
 z_Kc=PQ(2:K);
%
 zc=zeros(size(t));
 i=1; % start from interval 1: (t_1,t_2)
 for j=1:length(t)
     zc(j)=z_Kc(i);
     if t(j)>t_K(i+1)
        i=i+1;
     end
 end
%
%
 XL=levant_step(h_K,y_K,L2num);
%
% eps_hgo=0.01; % Fig. 1a
 eps_hgo=0.001; % Fig. 1b and 2a
% eps_hgo=0.02; % Fig. 2b
% eps_hgo=0.05; % Fig. 3
% eps_hgo=0.1;
 Xhgo=hgo_step(h_K,y_K,L1num,eps_hgo);
%
 figure,
 plot(t,dx,'r','LineWidth',2),hold on,
 plot(t_K,dx_K,'r.',...
 t_K,XL,'c.--',t_K,Xhgo,'g.--',...
    t,z,'k',t_K,z_K,'k.',... 
    t,zc,'b',t_K(1:end-1),z_Kc,'b.','LineWidth',1), grid,
legend_strings = { ...
    '$z(t)=dx/dt$ (analytical)', ...
    '', ...
    'Levant diff. (super-twisting)', ...
    ['HGO with $\varepsilon=' num2str(eps_hgo) '$'], ...
    ['quadratic spline with $\lambda=' num2str(lambda) '$'], ...
    '', ...
    ['zero-order spline with $\lambda=' num2str(lambda) '$'], ...
    '', ...
};
legend(legend_strings, ...
    'Location', 'south', ...
    'Interpreter', 'latex');
  xlabel('$t$','Interpreter','latex'), ylabel('$z(t)$ and $\hat z(t)$','Interpreter','latex'),
title( ...
  ['Approximation of derivatives with average sampling $h=' num2str(h) ...
   '$ and noise variance $\sigma=' num2str(sigma) '$'], ...
  'Interpreter', 'latex')

 disp(' '), disp('Variances:'), disp(' ')
% 
t_tr=floor(length(XL)/3);
disp(' - Levant'), disp(sqrt(mean((XL(t_tr:end)-dx_K(t_tr:end)).^2))),
disp(' - HGO'), disp(sqrt(mean((Xhgo(t_tr:end)-dx_K(t_tr:end)).^2))),
disp(' - zero-order Spline'), disp(sqrt(mean((z_Kc(t_tr-1:end)-dx_K(t_tr:end)).^2))),
disp(' - Quadratic Spline'), disp(sqrt(mean((z_K(t_tr:end)-dx_K(t_tr:end)).^2))),

%return
%
 figure,
 subplot(221),
 plot(t,dx,'r',t_K,dx_K,'r.','LineWidth',3), hold on,
 plot(t,z,'k',t_K,z_K,'k.','LineWidth',0.5), grid,
legend( ...
    '$z(t)=dx/dt$ (analytical)', ...
    '', ... % '$dx/dt$ at $t_k$', ...
    ['quadratic spline with $\lambda=' num2str(lambda) '$'], ...
    '', ... % 'quadratic spline at $t_k$', ...
    'Location', 'south', 'Interpreter', 'latex');
  xlabel('$t$','Interpreter','latex'), ylabel('$z(t)$ and $\hat z(t)$','Interpreter','latex'),
    title(['Approximation of derivatives with average sampling $h=' num2str(h) '$'], 'Interpreter', 'latex')
%
 subplot(222),
 plot(t,dx,'r','LineWidth',3),hold on,
 plot(t_K,dx_K,'r.',...
    t,zc,'b',t_K(1:end-1),z_Kc,'b.','LineWidth',1), grid,
legend( ...
    '$z(t)=dx/dt$ (analytical)', ...
    '', ... % '$dx/dt$ at $t_k$', ...
    ['zero-order spline with $\lambda=' num2str(lambda) '$'], ...
    '', ... % 'quadratic spline at $t_k$', ...
    'Location', 'south', 'Interpreter', 'latex');
  xlabel('$t$','Interpreter','latex'), ylabel('$z(t)$ and $\hat z(t)$','Interpreter','latex'),
title(['and noise variance $\sigma=' num2str(sigma) '$'], 'Interpreter', 'latex')
%
 subplot(223),
 plot(t,dx,'r','LineWidth',3),hold on,
 plot(t_K,dx_K,'r.',...
 t_K,XL,'c.--','LineWidth',1), grid,
 legend('$z(t)=dx/dt$ (analytical)','',...%'$dx/dt$ at $t_k$',...
    'Levant diff. (super-twisting)','Location','south','Interpreter','latex');
  xlabel('$t$','Interpreter','latex'), ylabel('$z(t)$ and $\hat z(t)$','Interpreter','latex'),
    title(['Approximation of derivatives with average sampling $h=' num2str(h) '$'], 'Interpreter', 'latex')

%
 subplot(224),
 plot(t,dx,'r','LineWidth',3),hold on,
 plot(t_K,dx_K,'r.',...
t_K,Xhgo,'g.--','LineWidth',1), grid,
 legend('$z(t)=dx/dt$ (analytical)','',...%'$dx/dt$ at $t_k$',...
  ['HGO with $\varepsilon=' num2str(eps_hgo) '$'],...
    'Location','south','Interpreter','latex');
  xlabel('$t$','Interpreter','latex'), ylabel('$z(t)$ and $\hat z(t)$','Interpreter','latex'),
title(['and noise variance $\sigma=' num2str(sigma) '$'], 'Interpreter', 'latex')
