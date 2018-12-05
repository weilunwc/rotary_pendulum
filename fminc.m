clear all
close all
clc

n=101;
dt = 0.01;
t = 0:dt:n*dt-dt;
x0 = [0;0;0;0];
xf = [0;pi;0;0];
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


xr = reshape(x_sol,[4,n]);
plot(t, u_sol)
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
plot(t,xr(2,:))
axis([0 max(t) -5 5])




T = 0:dt:10;
add_u = zeros(size(T,2),1);
add_u(1:size(u_sol,1))= u_sol;
ctr = timeseries(add_u,T);