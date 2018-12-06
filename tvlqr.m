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

[sol,fval,exitflag,output]  = fmincon(cost, initial_guess, [],[],[],[],lb,ub,@(x)constraint_z(x,n,dt))
x_sol = sol(1:4*n),
u_sol = sol(4*n+1:end);

xr = reshape(x_sol,[4,n]);

x_nom =xr';
u_nom = u_sol;

Q = blkdiag(10,10,1,1);
R = 1;
x_des = zeros(size(x_nom));
x_des(:,3) = linspace(-pi,0,size(x_nom,1));
x_des(:,4) = (pi/time(end))*ones(size(x_nom,1),1);
u_des = zeros(size(u_nom));



Qf = 100*eye(4);
S2T = Qf(:);
S1T = (-2*Qf*(x_des(end,:)-x_nom(end,:))');
S0T = (x_des(end,:)-x_nom(end,:))*Qf*(x_des(end,:)-x_nom(end,:))';
ST = [S2T;S1T;S0T];
[emit,Sflipped] = ode45(@(t,y)Sdynamics(t,y,Q,R,x_des,u_des,x_nom,u_nom,time),[time(end) 0],ST);
time_samples = flip(emit);
S_samples = flip(Sflipped); 

IC = [0;-pi;0;0];
[times,x_sol] =  ode45(@(t,x)optimal_dynamics(t,x,time,x_des,u_des,x_nom, u_nom,time_samples,S_samples),[0 time(end)],IC);
xr = interp1(times,x_sol,time)'
x_sim = xr;
if PLOT == 1
    plot(time, u_sol)
    xlabel('Step')
    ylabel('Control')
    figure()
    plot(xr(1,:),xr(3,:))
    xlabel('X-Position')
    ylabel('X-Velocity')
    figure()
    plot(xr(2,:),xr(4,:))
    xlabel('Angle[rad]')
    ylabel('Angular Velocity[rad/s]')
    figure()
    xlabel('Time')
    ylabel('pos')
    plot(times,x_sol(:,2))
    axis([0 max(times) -5 5])
end
