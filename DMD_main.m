clear; close all; clc;

%% Define time and space discretizations
xi = linspace(-10,10,400);
t = linspace(0,4*pi,200); 
dt = t(2) - t(1);
[Xgrid,T] = meshgrid(xi,t);

%% Create two spatio-temporal patterns
f1 = sech(Xgrid+3) .* (1*exp(1j*2.3*T));
f2 = (sech(Xgrid).*tanh(Xgrid)).*(2*exp(1j*2.8*T));

%% Combine signals and make data matrix
f = f1 + f2;
X = f.'; % Data Matrix

%% Create DMD data matrices and ranks
X1 = X(:, 1:end-1);
X2 = X(:, 2:end);
r = 2; % rank truncation

%% Compute DMD Solution
[Phi,omega,lambda,b,X_dmd] = DMD(X1,X2,r,dt);

%% Plotting
fig1 = figure(1);
fig1.Color = 'w';
fig1.Position = [300,300,700,300];

subplot(1,2,1); 
surfl(real(f)); 
shading interp; colormap(gray); view(-20,60);
xlabel('x');
ylabel('t');
title('Original dynamics')
set(gca, 'YTick', numel(t)/4 * (0:4)), 
set(gca, 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
set(gca, 'XTick', linspace(1,numel(xi),3)), 
set(gca, 'Xticklabel',{'-10', '0', '10'});

subplot(1,2,2); 
surfl(real(X_dmd')); 
shading interp; colormap(gray); view(-20,60);
xlabel('x');
ylabel('t');
title('DMD reconstruction')
set(gca, 'YTick', numel(t)/4 * (0:4)), 
set(gca, 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
set(gca, 'XTick', linspace(1,numel(xi),3)), 
set(gca, 'Xticklabel',{'-10', '0', '10'});
