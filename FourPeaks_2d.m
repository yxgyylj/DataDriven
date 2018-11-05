%%--------------------- Simulating Gray Scott model -------------------------
%
%   Model discription:
%       G-S model in 1d with 4 peaks (coral-like sturcture in 2d)
%
%	Written by Xige Yang
%   Ohio State University
%   
%   --  created on 04/22/2018
%       first revised on 04/25/2018
%-------------------------------------------------------------------------
clear; close all; clc;
% Feed and kill rates
rho=.055;
mu=.062;
% rho=.01;
% mu=1;

% Diffusion rates
DS = 1.;
DA = .1;

% Domain size and grid size
L = 1;
Nx = 256;    dx = 1/Nx;
Ny = 256;    dy = 1/Ny;
x = linspace(0,L-dx,Nx);
DS = DS * dx^2; DA = DA * dx^2; % we need to rescale diffusion rates

% 5,000 simulation seconds with 4 steps per simulated second
t0 = 0;
dt = .25;
T = 2500;
tspan = t0:dt:T;
Nt = length(t0+dt:dt:T);
[S, A] = ICs(Nx, 1);

% save data
save_step = 20;
save_Nt = 0:save_step:Nt;
save_tspan = floor(save_Nt*dt);
[xx, tt] = meshgrid(x,save_tspan);
N_save = floor(Nt/save_step)+1;
S_save = zeros(Nx, N_save);
A_save = zeros(Nx, N_save);
S_save(:,1) = S;
A_save(:,1) = A;

% Add plot vectors
normA = zeros(Nt+1,1);
DirechA = zeros(Nt+1,1);
normA(1) = norm(A);
DirechA(1) = norm(MyLaplacian1(A));

% Add a scaled-color image
hi = plot(1:Nx,A);
s=get(gcf, 'Position');
s(3)=400;
s(4)=300;
set(gcf, 'Position', s);
%ylim([0,1]);

% Text setup
ht = text(3,Nx-3,'Time = 0');
ht.Units = 'normalized';
ht.Position=[.35,1.05];
ht.Color = [.95 .2 .8];
ht.FontSize = 16;
drawnow

% Start simulation
tic
nframes = 1;
for j = 1:Nt
    % Atomatically used periodic boundary condition
    anew = S + (DS*MyLaplacian1(S)/dx^2 - S.*A.^2 + rho*(1-S))*dt;
    bnew = A + (DA*MyLaplacian1(A)/dx^2 + S.*A.^2 - (mu+rho)*A)*dt;
    S = anew;
    A = bnew;
    hi.YData = A;
    ht.String = ['Time = ' num2str(j*dt)];
    
    if(mod(j,save_step) == 0)
        nframes = nframes + 1;
        S_save(:,nframes) = S;
        A_save(:,nframes) = A;
        drawnow
    end
    
    normA(j+1) = norm(A);
    DirechA(j+1) = norm(MyLaplacian1(A));
end

delta = toc;
disp([num2str(nframes) ' frames in ' num2str(delta) ' seconds']);

% save data
% save('Data/solutions_GS_1d_onepeak.mat','S_save','A_save','tt','xx'...
%     ,'save_step','dt', 'x', 'save_tspan');

axes('Position',[0 0 1 1])
axis off

% plot the outcome
fig2 = figure;
subplot(2,1,1)
plot(tspan(80:end), normA(80:end));
xlabel('time', 'FontSize', 24);
ylabel('\int_{\Omega} A^2(x,y) dxdy', 'FontSize', 24)
title('A concentration', 'FontSize', 28);
subplot(2,1,2)
plot(tspan(80:end), DirechA(80:end));
xlabel('time', 'FontSize', 24);
ylabel('C \int_{\Omega} (\nabla A(x,y))^2 dxdy', 'FontSize', 24)
title('Rescaled Direchlet energy of A', 'FontSize', 28);
% saveas(fig2,'B_plot.png');

figure('Position',[300,300,700,300]);
subplot(1,2,1)
surf(xx,tt,S_save');
xlabel('x');
ylabel('t');
ylim([t0,T]);
title('S');
shading interp;

subplot(1,2,2)
surf(xx,tt,A_save');
xlabel('x');
ylabel('t');
ylim([t0,T]);
title('A');
shading interp;

% Add a scaled-color image
%hi = image(B);
%hi.CDataMapping = 'scaled';

function [S, A] = ICs(m,n)
  % Initialize A to one
  S = ones(m,n);
  % Initialize B to zero which a clump of ones
  A = zeros(m,n);
  oneThird = floor(m/3);
  A(oneThird+1:2*oneThird) = 1;
end

function out = MyLaplacian1(in)
  out = - 2*in ...
      + (circshift(in,[ 1, 0]) + circshift(in,[-1, 0]));
end