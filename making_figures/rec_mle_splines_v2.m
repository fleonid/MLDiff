% FILE: ml_rec_mle_splines_v2.m
% DESCRIPTION:
%   Runs the full pipeline for maximum-likelihood-based recursive quadratic
%   and zero-order spline differentiation, including data generation, spline fitting, 
%   and comparison with other methods.
%
% CONTEXT:
%   Based on the method described in:
%   "Numerical Differentiation under Coarse Non-uniform Sampling and Gaussian Noise"
%   by K. E. Avrachenkov and L. B. Freidovich, 2025.
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
%   2025-07-28
%
% direct='true';
 direct='false';
%
 diary rec_mle_splines_v2.txt
 disp(' '), disp(datetime('now')), %disp(' '),
%
%switch direct
%    case 'true'
         clearvars
         close all
%         direct='true';
          direct='false';
%    otherwise
         disp('recursive computations')
%end
%
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultAxesLineWidth', 1.5); 
set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultTextFontSize', 14);
set(0, 'DefaultFigureWindowState', 'maximized');
%
 fig_case = 'Fig2a';  % options: 'Fig1a', 'Fig1b', 'Fig2a', 'Fig2b', 'Fig3'
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
% h=0.0002; % very small sampling / long computation time
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
%switch direct
%    case 'true'
 [K,t_K,h_K,x_K,dx_K,y_K,L0num,L1num,L2num,t,x,dx]=make_data(h,dh,t_initial,t_final,sigma);
%    otherwise
%         disp('using previously generated data'),
%end
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
    lambda0=0.01;  % special case of very noisy data
end
%
 disp(' '), fprintf('lambda1 = %g', lambda1), disp(' ') % for quadratic splines
 disp(' '), fprintf('lambda0 = %g', lambda0), disp(' ') % for 0-order splines
%
 K0=5;
%
%  QUADRATIC SPLINES
%
 tic,
%
%  Initialization for recursive computations
%
 hK=h_K(1:K0-1);
 yK=y_K(1:K0);
% 
 [Q,C]=quadratic_spline_step_QC(hK); 
 A=C'*C+lambda1*Q;
 invA=inv(A);
 PQ=invA*C'*yK; 
 x0=PQ(1); 
 pK=PQ(2:K0);
 zK=z_from_p(pK,hK);
%
 p_K=pK;
 z_K=zK;
%
%  ZERO ORDER SPLINES
%
 [Qc,Cc]=zero_order_spline_step_QC(hK); 
 Ac=Cc'*Cc+lambda0*Qc;
 invAc=inv(Ac);
 PQ=invAc*Cc'*yK;
 Hc=[1;hK];
 Bc=yK(1);
 for i=1:K0-1
   Bc=[Bc;0]; Bc=Bc+Hc(1:i+1)*yK(i+1);
 end
%
 x0c=PQ(1); 
 zKc=PQ(2:K0);
%
 z_Kc=zKc;
%
 for Ki=K0+1:K
%
 hK=h_K(1:Ki-1);
 yK=y_K(1:Ki);
% 
%  QUADRATIC SPLINES
%
switch direct
case 'true'
%
% not recursive
% 
 [Q,C]=quadratic_spline_step_QC(hK); 
 PQ=(C'*C+lambda1*Q)\C'*yK; 
 x0=PQ(1); 
 pK=PQ(2:Ki);
otherwise
%
% recursive
%
 [dx0,dP,pKp,Q,C,invA]=update_quadratic_spline(x0,pK,hK,yK(end-1),yK(end),lambda1,Q,C,invA);
 x0=x0+dx0;
 pK=[pK+dP;pKp];
%
end
%
 zK=z_from_p(pK,hK);
 p_K=[p_K;pK(end)];
 z_K=[z_K;zK(end)];
%
%  ZERO ORDER SPLINES
%
switch direct
case 'true'
%
% not recursive
% 
 [Qc,Cc]=zero_order_spline_step_QC(hK); 
 PQ=(Cc'*Cc+lambda0*Qc)\Cc'*yK; 
 x0c=PQ(1); 
 zKc=PQ(2:Ki);
%
otherwise
%
% recursive
%
 Hc=[Hc;hK(end)];
 Bc=[Bc;0]; 
 Bc=Bc+Hc*yK(end);
%
 invAs=[invAc,zeros(Ki-1,1);zeros(1,Ki-1),1/(lambda0*hK(end))];
 invAc=invAs-invAs*Hc*transpose(Hc)*invAs/(1+transpose(Hc)*invAs*Hc);
 zKc=invAc*Bc;
%
end
%
 z_Kc=[z_Kc;zKc(end)];
%
% t_K(Ki),
%
 end
%
%  Interpolations
%
 z=interp_z_from_z_K(t,t_K,h_K,z_K,p_K);
%
 zc=zeros(size(t));
 i=1; % start from interval 1: (t_1,t_2)
 for j=1:length(t)
     zc(j)=z_Kc(i);
     if t(j)>t_K(i+1)
        i=i+1;
     end
 end
 disp(' '), toc, disp(' '), 
%
%  Super twisting differentiator
%
 XL=levant_step(h_K,y_K,L2num);
%
% High-gain observer
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
%  Plotting figures
%
%  figure,
%  plot(t,dx,'r','LineWidth',3),hold on,
%  plot(t_K,dx_K,'r.',...
%  t_K,XL,'c.--',t_K,Xhgo,'g.--',...
%     t,z,'k',t_K,z_K,'k.',... 
%     t,zc,'b',t_K(1:end-1),z_Kc,'b.','LineWidth',2), grid,
% legend_strings = { ...
%     '$z(t)=dx/dt$ (analytical)', ...
%     '', ...
%     'Levant diff. (super-twisting)', ...
%     ['HGO with $\varepsilon=' num2str(eps_hgo) '$'], ...
%     ['quadratic spline with $\lambda=' num2str(lambda1) '$'], ...
%     '', ...
%     ['zero-order spline with $\lambda=' num2str(lambda0) '$'], ...
%     '', ...
% };
% legend(legend_strings, ...
%     'Location', 'south', ...
%     'Interpreter', 'latex');
%   xlabel('$t$','Interpreter','latex'), ylabel('$z(t)$ and $\hat z(t)$','Interpreter','latex'),
% title( ...
%   ['Approximation of derivatives with average sampling $h=' num2str(h) ...
%    '$ and noise variance $\sigma=' num2str(sigma) '$'], ...
%   'Interpreter', 'latex')

 disp(' '), disp('RMSEs:'), disp(' ')
% 
%  Computing Root mean squuare errors after transients (last two-thirds)
%
t_tr=floor(length(XL)/3);
disp(' - Super Twisting'), disp(sqrt(mean((XL(t_tr:end)-dx_K(t_tr:end)).^2))),
disp(' - HGO'), disp(sqrt(mean((Xhgo(t_tr:end)-dx_K(t_tr:end)).^2))),
disp(' - Zero-Order Spline'), disp(sqrt(mean((z_Kc(t_tr-1:end)-dx_K(t_tr:end)).^2))),
disp(' - Quadratic Spline'), disp(sqrt(mean((z_K(t_tr:end)-dx_K(t_tr:end)).^2))),
%
 diary off
%
%  Plotting figures
%
 figure, axis equal
 subplot(221),
 plot(t,dx,'r',t_K,dx_K,'r.','LineWidth',3), hold on,
 plot(t,z,'k',t_K,z_K,'k.','LineWidth',2), grid,
legend( ...
    '$z(t)=dx/dt$ (analytical)', ...
    '', ... % '$dx/dt$ at $t_k$', ...
    ['quadratic spline with $\lambda=' num2str(lambda1) '$'], ...
    '', ... % 'quadratic spline at $t_k$', ...
    'Location', 'south', 'Interpreter', 'latex');
  xlabel('$t$','Interpreter','latex'), ylabel('$z(t)$ and $\hat z(t)$','Interpreter','latex'),
    title(['Approximation of derivatives with average sampling $h=' num2str(h) '$'], 'Interpreter', 'latex')
    ylim([-2,2]);
%
 subplot(222),
 plot(t,dx,'r','LineWidth',3),hold on,
 plot(t_K,dx_K,'r.',...
    t,zc,'b',t_K(1:end-1),z_Kc,'b.','LineWidth',2), grid,
legend( ...
    '$z(t)=dx/dt$ (analytical)', ...
    '', ... % '$dx/dt$ at $t_k$', ...
    ['zero-order spline with $\lambda=' num2str(lambda0) '$'], ...
    '', ... % 'quadratic spline at $t_k$', ...
    'Location', 'south', 'Interpreter', 'latex');
  xlabel('$t$','Interpreter','latex'), ylabel('$z(t)$ and $\hat z(t)$','Interpreter','latex'),
 title(['and noise variance $\sigma=' num2str(sigma) '$'], 'Interpreter', 'latex')
 ylim([-2,2]);
%
 subplot(223),
 plot(t,dx,'r','LineWidth',3),hold on,
 plot(t_K,dx_K,'r.',...
 t_K,XL,'c.--','LineWidth',2), grid,
 legend('$z(t)=dx/dt$ (analytical)','',...%'$dx/dt$ at $t_k$',...
    'Levant diff. (super-twisting)','Location','south','Interpreter','latex');
  xlabel('$t$','Interpreter','latex'), ylabel('$z(t)$ and $\hat z(t)$','Interpreter','latex'),
    title(['Approximation of derivatives with average sampling $h=' num2str(h) '$'], 'Interpreter', 'latex')

%
 subplot(224),
 plot(t,dx,'r','LineWidth',3),hold on,
 plot(t_K,dx_K,'r.',...
t_K,Xhgo,'g.--','LineWidth',2), grid,
 legend('$z(t)=dx/dt$ (analytical)','',...%'$dx/dt$ at $t_k$',...
  ['HGO with $\varepsilon=' num2str(eps_hgo) '$'],...
    'Location','south','Interpreter','latex');
  xlabel('$t$','Interpreter','latex'), ylabel('$z(t)$ and $\hat z(t)$','Interpreter','latex'),
title(['and noise variance $\sigma=' num2str(sigma) '$'], 'Interpreter', 'latex')
