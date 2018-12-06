clear all
close all
clc

PLOT = 1

n=101;
dt = 0.01;
time = 0:dt:n*dt-dt;
x0 = [0;-pi;0;0];
xf = [0;0;0;0];
initial_guess = [x0;ones((n-2)*4,1);xf;ones(n,1)]; %[all states except, all actions]

Z = zeros(n*4,1);
Z(1:4:end) = -3*pi/4;
Z(2:4:end) = -2*pi;
Z(3:4:end) = -Inf;
Z(4:4:end) = -Inf;
Z(1:4) = x0;
Z((n-1)*4+1:n*4) = xf;
W = zeros(n*4,1);
W(1:4:end) = 3*pi/4;
W(2:4:end) = 2*pi;
W(3:4:end) = Inf;
W(4:4:end) = Inf;
W(1:4) = x0;
W((n-1)*4+1:n*4) = xf;

lb = [Z;-Inf*ones(n,1)];
ub = [W; Inf*ones(n,1)];
cost = @(x) x(4*n+1:end)'*x(4*n+1:end); 

[sol,fval,exitflag,output]  = fmincon(cost, initial_guess, [],[],[],[],lb,ub,@(x)constraint(x,n,dt))
x_sol = sol(1:4*n),
u_sol = sol(4*n+1:end);


x_nom = reshape(x_sol,[4,n]);
u_nom = u_sol;

Q = blkdiag(10,10,1,1);
R = 1;

Qf = 100*eye(4);
% S2T = Qf(:);
% S1T = (-2*Qf*(x_nom(end,:)'));
% S0T = (x_nom(end,:))*Qf*(x_nom(end,:))';
% ST = [S2T;S1T;S0T];
% [emit,Sflipped] = ode45(@(t,y)Snomdynamics(t,y,Q,R,x_nom,u_nom,time),[time(end) 0],ST);
% time_samples = flip(emit);
% S_samples = flip(Sflipped); 



IC = [0;-pi;0;0];
x(:,1) = IC;

for i = 1:length(time)
[A,B] = linearizedfurata_z(x(:,i),0);
[tflipped,Sflipped] = ode45(@(t,S)SABdynamics(t,S,A,B,Q,R),[time(end) 0],Qf);
t = flip(tflipped);
S = flip(Sflipped);
%[t,x] = ode45(@(t,x)dynamics(t,x,A,B,R,S,xd,time),[0 T],x0);

S_at_t = interp1(t,S,time(1,i));
Su = [S_at_t(1:4)' S_at_t(5:8)' S_at_t(9:12)' S_at_t(13:16)']; 
xn = x_nom(:,i);
un = u_nom(i);

u(i) = un-inv(R)*B'*Su*(x(:,i)-xn);
fnom = furataDynamics_z(xn,un);
xdot = fnom + A*(x(:,i)-xn)+B*(u(i)-un);
x(:,i+1) = x(:,i) + xdot*dt;
end

u_sol = u;
x_sol = x(:,1:length(time));

%[times,x_sol] =  ode45(@(t,x)optimal_dynamics_basic(t,x,time,x_nom, u_nom,time_samples,S_samples),[0 time(end)],IC);


if PLOT == 1
    plot(time, u_sol)
    xlabel('Step')
    ylabel('Control')
    figure()
    plot(x_sol(1,:),x_sol(3,:))
    xlabel('X-Position')
    ylabel('X-Velocity')
    figure()
    plot(x_sol(2,:),x_sol(4,:))
    xlabel('Angle[rad]')
    ylabel('Angular Velocity[rad/s]')
    figure()
    xlabel('Time')
    ylabel('pos')
    plot(time,x_sol(2,:))
    axis([0 max(time) -5 5])
end

x_sim = x_sol;