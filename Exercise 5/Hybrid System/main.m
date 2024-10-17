%% Simulation of bouncing ball

clc; clear; close all;

% initial conditions
x0 = [1 0];

% simulation horizon
TSPAN=[0 10];
JSPAN = [0 20];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',0.1);

% simulate
[t,j,x] = HyEQsolver( @f_ex1,@g_ex1,@C_ex1,@D_ex1,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% plot solution
figure(1) % position
clf
subplot(2,1,1), plotHarc(t,j,x(:,1));
grid on
xlabel('Time (s)')
ylabel('$x_1$','Interpreter','Latex','Fontsize',12)
title('position')
subplot(2,1,2), plotHarc(t,j,x(:,2));
grid on
xlabel('Time (s)')
ylabel('$x_2$','Interpreter','Latex','Fontsize',12)
title('velocity')

% plot phase plane
figure(2) % position
clf
plotHarcColor(x(:,1),j,x(:,2),t);
xlabel('$x_1$','Interpreter','Latex','Fontsize',12)
ylabel('$x_2$','Interpreter','Latex','Fontsize',12)
title('phase portrait')
grid on

% plot hybrid arc
figure(3)
plotHybridArc(t,j,x)
axis([-inf inf -inf inf 0 1.5])
grid on
xlabel('$j$','Interpreter','Latex','Fontsize',12)
ylabel('$t$','Interpreter','Latex','Fontsize',12)
zlabel('$x_1$','Interpreter','Latex','Fontsize',12)

