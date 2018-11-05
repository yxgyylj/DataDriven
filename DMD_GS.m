clear; close all; clc;
set(0,'Defaultaxesfontsize',16)

%% Define time and space discretizations
load('Data/solutions_GS_1d.mat')


%% Create DMD data matrices and ranks
S_D1 = S_save(:, 1:end-1);
S_D2 = S_save(:, 2:end);
A_D1 = A_save(:, 1:end-1);
A_D2 = A_save(:, 2:end);
dt = dt*save_step;
r = 128; % rank truncation
r1 = 64;
r2 = 32;

%% Compute DMD Solution
[Phi_S,omega_S,lambda_S,b_S,S_dmd] = DMD(S_D1,S_D2,r,save_step);
[Phi_A,omega_A,lambda_A,b_A,A_dmd] = DMD(A_D1,A_D2,r,save_step);
[Phi_S_1,omega_S_1,lambda_S_1,b_S_1,S_dmd_1] = DMD(S_D1,S_D2,r1,save_step);
[Phi_A_1,omega_A_1,lambda_A_1,b_A_1,A_dmd_1] = DMD(A_D1,A_D2,r1,save_step);
[Phi_S_2,omega_S_2,lambda_S_2,b_S_2,S_dmd_2] = DMD(S_D1,S_D2,r2,save_step);
[Phi_A_2,omega_A_2,lambda_A_2,b_A_2,A_dmd_2] = DMD(A_D1,A_D2,r2,save_step);

%% Plotting
fig1 = figure(1);
fig1.Color = 'w';
fig1.Position = [300,300,500,500];

subplot(2,2,1)
imagesc(save_tspan,x,S_save);
xlabel('t');
ylabel('x');
title('Original S');

subplot(2,2,2)
imagesc(save_tspan,x,A_save);
xlabel('t');
ylabel('x');
title('Original A');

subplot(2,2,3); 
imagesc(save_tspan,x,real(S_dmd));
xlabel('t');
ylabel('x');
title(sprintf('DMD at rank = %d',r));

subplot(2,2,4);
imagesc(save_tspan,x,real(A_dmd));
xlabel('t');
ylabel('x');
title(sprintf('DMD at rank = %d',r));

fig2 = figure(2);
fig2.Position = [300,300,500,500];

subplot(2,2,1)
surf(tt(2:end,:),xx(2:end,:),S_save(:,2:end)');
shading interp; view(-20,60);xlim([save_tspan(2),save_tspan(end)]);
xlabel('t');
ylabel('x');
title('Original S');

subplot(2,2,2)
surf(tt(2:end,:),xx(2:end,:),A_save(:,2:end)');
shading interp; view(-20,60);xlim([save_tspan(2),save_tspan(end)]);
xlabel('t');
ylabel('x');
title('Original A');

subplot(2,2,3); 
surf(tt(2:end,:),xx(2:end,:),real(S_dmd)');
shading interp; view(-20,60);xlim([save_tspan(2),save_tspan(end)]);
xlabel('t');
ylabel('x');
title(sprintf('DMD at rank = %d',r));

subplot(2,2,4);
surf(tt(2:end,:),xx(2:end,:),real(A_dmd)');
shading interp; view(-20,60);xlim([save_tspan(2),save_tspan(end)]);
xlabel('t');
ylabel('x');
title(sprintf('DMD at rank = %d',r));


%% Play movie
hi = figure('Position',[300,300,300,200]);
count = 1;
plt1 = plot(x,real(A_dmd(:,1)),'linewidth',2);
hold on;
plt2 = plot(x,real(A_dmd_1(:,1)),'linewidth',2);
plt3 = plot(x,real(A_dmd_2(:,1)),'linewidth',2);
plt4 = plot(x,A_save(:,1),'linewidth',2);
ylim([0,1]);
legend(sprintf('rank = %d',r),sprintf('rank = %d',r1),...
    sprintf('rank = %d',r2),'Original S');
tit = title(sprintf('Time = %.3f', save_tspan(1)));

%% Plot error
gif_ts = 10:15:length(save_tspan);
errors = zeros(3,1+length(gif_ts));
errors(:,1) = [norm(real(A_dmd(:,1))-A_save(:,1));...
                norm(real(A_dmd_1(:,1))-A_save(:,1));...
                norm(real(A_dmd_2(:,1))-A_save(:,1))];
saveas(hi,sprintf('Imgs/GS_%03d.png',count));
% error

for j = gif_ts
    count = count + 1;
    plt1.YData = real(A_dmd(:,j));
    plt2.YData = real(A_dmd_1(:,j));
    plt3.YData = real(A_dmd_2(:,j));
    plt4.YData = A_save(:,j);
    tit.String = sprintf('Time = %.3f', save_tspan(j));
    errors(:,1) = [norm(real(A_dmd(:,j))-A_save(:,j));...
                norm(real(A_dmd_1(:,j))-A_save(:,j));...
                norm(real(A_dmd_2(:,j))-A_save(:,j))];
    saveas(hi,sprintf('Imgs/GS_%03d.png',count));
    drawnow;
end
figure('Position',[300,300,400,300]);
plot(0:length(gif_ts),errors);
legend(sprintf('rank = %d',r),sprintf('rank = %d',r1),sprintf('rank = %d',r2));
