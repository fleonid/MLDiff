% FILE: ml_quadratic_spline_v16.m
% DESCRIPTION:
%   Runs the full pipeline for maximum-likelihood-based quadratic spline differentiation, including data generation, spline fitting, and comparison with other methods.
%
% CONTEXT:
%   Based on the method described in:
%   "Numerical Differentiation under Coarse Non-uniform Sampling and Gaussian Noise"
%   by K. E. Avrachenkov and L. B. Freidovich, submitted to IEEE TAC, 2025.
%
% INPUTS:
%   None (parameters are set within the script)
%
% OUTPUTS:
%   Plots and printed performance metrics
%
% AUTHOR:
%   Leonid B. Freidovich (leonid.freidovich@umu.se)
%
% LAST UPDATED:
%   2025-07-22
 clearvars
 close all
%
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultAxesLineWidth', 1.5); 
set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultTextFontSize', 14);
%
% diary ml_quadratic_v16.txt
 disp(' '), disp(datetime('now')),
%
 fig_case = 'Fig3';  % options: 'Fig1a', 'Fig1b', 'Fig2a', 'Fig2b', 'Fig3'
%
% standard deviation of Gaussian noise
%
switch fig_case
    case {'Fig1a', 'Fig1b', 'Fig3'}
        sigma = 1e-4;
    case 'Fig2a'
        sigma = 1e-7;
    case 'Fig2b'
        sigma = 1e-2;
    otherwise
        sigma = 1e-3;  % default
end
 disp(' '), fprintf('standard deviation of Gaussian noise sigma = %f', sigma), disp(' '),
%
%    average sampling step
%
switch fig_case
    case 'Fig1a'
        h       = 0.01;
    case {'Fig1b', 'Fig2a','Fig2b'}
        h       = 0.001;
    case 'Fig3'
        h       = 0.05;
    otherwise
        h       = 0.075;
end
disp(' '), fprintf('average sampling  h = %f', h), disp(' '),
%
 t_initial=0;
 t_final=1.95;
%
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
 lambda1=1e-4; % empirically found
 lambda0=1e-4;
% lambda=2.2e-7; % theoretically optimal
%
if  strcmp(fig_case, 'Fig2b')
    lambda0=0.1;
end
%
 disp(' '), fprintf('lambda1 = %g', lambda1), disp(' ')
 disp(' '), fprintf('lambda0 = %g', lambda0), disp(' ')
%
%  QUADRATIC SPLINES
%
% [x0,p_K]=quadratic_spline_step(h_K,y_K,lambda1);
 [Q,C]=quadratic_spline_step_QC(h_K); 
 PQ=(C'*C+lambda1*Q)\C'*y_K; 
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
 PQ=(C'*C+lambda0*Q)\C'*y_K; 
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
switch fig_case
    case 'Fig1a'
        eps_hgo = 0.01;
    case {'Fig1b', 'Fig2a'}
        eps_hgo = 0.001;
    case 'Fig2b'
        eps_hgo = 0.02;
    case 'Fig3'
        eps_hgo = 0.05;
    otherwise
        eps_hgo = 0.1;  % fallback value
end
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
    ['quadratic spline with $\lambda=' num2str(lambda1) '$'], ...
    '', ...
    ['zero-order spline with $\lambda=' num2str(lambda0) '$'], ...
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
 plot(t,z,'k',t_K,z_K,'k.','LineWidth',1), grid,
legend( ...
    '$z(t)=dx/dt$ (analytical)', ...
    '', ... % '$dx/dt$ at $t_k$', ...
    ['quadratic spline with $\lambda=' num2str(lambda1) '$'], ...
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
    ['zero-order spline with $\lambda=' num2str(lambda0) '$'], ...
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
