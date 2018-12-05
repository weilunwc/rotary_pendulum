% 16748 Underactuated Robotics Final Project

clear;
close all;
clc;

global m l g b
m = 1;
l = 1;
g = 1;
b = 0.1;
% Specify the target
x_target = [pi;0];
tol = 1e-4;
%%  Shooting method

% Use shooting method to solve swing up problem
torque_limit = 1;
upper_bound = torque_limit;
lower_bound = -torque_limit;
T = 10;
dt = 0.025;

tic
[u_traj, x_traj] = shooting_method(upper_bound, lower_bound, dt, T);
toc
x0 = [0;0];
x_traj = [x0 x_traj];
% Plot the trajectory
sita = x_traj(1,:);
sita_dot = x_traj(2,:);
subplot(1,2,1);
plot(sita, sita_dot,'LineWidth',3);
xlabel('\theta');
ylabel('$\dot{\theta}$','Interpreter','latex');
title('phase portrait');
set(gca,'fontsize',20);

subplot(1,2,2);
n = T/dt;
plot([1:n]*T/n,u_traj(1:n),'LineWidth',3);
title('control history');
xlabel('u');
set(gca,'fontsize',20);

%% animation
animation(sita);


%% RG-RRT implementation
%% 
tic
sita_range = [-4 4];
sita_dot_range = [-2 2];
start = [0;0];
goal = [-pi;0];
goal2 = [pi;0];
goal_bias = 0.1;

% Initialize tree
id = 1;
q_init.coord = [0;0];
q_init.parent = [];
q_init.id = 1;
q_init.control = [];
q_init.parent_id = 0;
[tmax, ymax] = ode45(@(t,y)pendulum_dynamics(t,y,torque_limit),[0,0.1], start);
[tmin, ymin] = ode45(@(t,y)pendulum_dynamics(t,y,-torque_limit),[0,0.1], start);
q_init.reachable_set = zeros(2,2);
q_init.reachable_set(:,1) = ymin(end,:)';
q_init.reachable_set(:,2) = ymax(end,:)';
id = id + 1;

iterations = 1000;
Tree.points = repmat(q_init,iterations+1,1);


for i = 1:iterations
    while true
        % goal bias
        is_goal_bias = false;
        if rand() <= goal_bias
            is_goal_bias = true;
            random_coord = goal;
        else
            sita = rand() * (sita_range(2) - sita_range(1)) + sita_range(1);
            sita_dot = rand() * (sita_dot_range(2) - sita_dot_range(1)) + sita_dot_range(1);
            random_coord = [sita;sita_dot];
        end
        % find nearest point in tree
        dist = inf;
        n_points = max(size(Tree.points));
        nearest_coord = [];
        for j = 1:n_points
            candidate = Tree.points(j).coord;
            candidate_id = Tree.points(j).id;
            if norm(random_coord - candidate, 2) < dist
                dist = norm(random_coord - candidate, 2);
                nearest_coord = candidate;  
                nearest_id = candidate_id;
            end
        end
        
        if is_goal_bias == true
            
            % Since this is goal biased, don't care if it expands
            % efficiently
            break;
        end
        
        % Check if the nearest point efficeiently expands the workspace
        nearest_point = Tree.points(nearest_id);
        
        if dist > norm(random_coord - nearest_point.reachable_set(:,1))
            break
        end
        if dist > norm(random_coord - nearest_point.reachable_set(:,2))
            break
        end      

    end
    % reach to the sampled point from the nearest coord
    dist = inf;
    for u = -torque_limit:torque_limit:torque_limit
        [ts, ys] = ode45(@(t,y)pendulum_dynamics(t,y,u),[0,0.1], nearest_coord);
        if norm(ys(end,:) - random_coord, 2) < dist
            dist = norm(ys(end,:) - random_coord,2);
            new_coord = ys(end,:)';
            best_u = u;
        end
    end
    
    % check if this point is new
    %{
    repeated = false;
    for j = 1:n_points
        candidate = Tree.points(j).coord;
        if new_coord(1) == candidate(1) && new_coord(2) == candidate(2)
            repeated = true;
            break;
        end
    end
    
    if repeated
        fprintf("repeated\n");
        continue;
    end
    %}
    q.coord = new_coord;
    q.parent = nearest_coord;
    q.parent_id = nearest_id;
    q.id = id;
    q.control = best_u;
    % Find reachable set
    [tmax, ymax] = ode45(@(t,y)pendulum_dynamics(t,y,torque_limit),[0,0.1], new_coord);
    [tmin, ymin] = ode45(@(t,y)pendulum_dynamics(t,y,-torque_limit),[0,0.1], new_coord);
    q.reachable_set = zeros(2,2);
    q.reachable_set(:,1) = ymin(end,:)';
    q.reachable_set(:,2) = ymax(end,:)';
    Tree.points(i+1) = q;
    
    % Check if reached goal
    if norm(new_coord - goal,2) < 0.1
        fprintf("reached goal at %d\n",i);
        break;
    end
end

% plot the tree
figure();
hold on;
for i = 1:n_points
    q = Tree.points(i);
    if size(q.parent) ~= 0
        plot([q.coord(1) q.parent(1)],[q.coord(2) q.parent(2)],'b');
    end
end

% plot the goal
scatter([-pi pi],[0 0],'black');
axis([sita_range sita_dot_range]);
title('RRT')
toc
%% 

% Track the path
% First find the goal
found_goal = false;
n_points = max(size(Tree.points));
for i = 1:n_points
    q = Tree.points(i);
    if norm(q.coord - goal) < 0.1
        goal_id = i;
        found_goal = true;
        break;
    end
end


% retrack the tree
if found_goal
    q = Tree.points(goal_id);
    path = [];
    while norm(q.coord-start,2) > 0.1 
        path = [path q.coord];
        parent_id = q.parent_id;
        q = Tree.points(parent_id);
    end

    path = [path q.coord];
    path = flip(path,2);

    % plot the path
    n_path = size(path,2);
    plot(path(1,:),path(2,:),'r','LineWidth',3);
end




function [u_traj, x_traj] = shooting_method(upper_bound, lower_bound, dt, T)
% Solve for the trajcetory optimization problem given control contraints

L = T/dt;
u0 = zeros(L,1);

% cost function
fun = @(u)u'*u; 
A = [];
b = [];
Aeq = [];
beq = [];


ub = upper_bound * ones(L,1);
lb = lower_bound * ones(L,1);

options =optimoptions(@fmincon,'TolFun',0.00000001,'MaxIter',10000,'MaxFunEvals',...
    100000,'Display','final','DiffMinChange', 0.001,'Algorithm','sqp');

nonlcon = @dynamics;


% Use shooting method to solve trajectory optimization
u_traj = fmincon(fun,u0,A,b,Aeq,beq,lb,ub,nonlcon,options);

% Forward pass the solved control trajectory into the real dynamics 
x_traj = forward_pass(u_traj);


end

function [c, ceq] = dynamics(u)
% The nonlinear dynamics constraint incorparated into the optimization
% problem
c = [];

x0 = [0;0];
x = x0;
global m l g b
dt = 0.025;
n = size(u);
for i = 1:n
    dx = zeros(2,1);
    dx(1) = x(2);
    dx(2) = (1/(m*l^2)) * (u(i) - b*x(2) - m*g*l*sin(x(1)));
    x(1) = x(1) + dt * dx(1);
    x(2) = x(2) + dt * dx(2);
end

ceq = [x(1) - pi; x(2) - 0];

end


function x_traj = forward_pass(u)
% Forward pass
global m l g b
x0 = [0;0];
x = x0;
dt = 0.025;
n = max(size(u));
x_traj = zeros(2,n);
for i = 1:n
    dx = zeros(2,1);
    dx(1) = x(2);
    dx(2) = (1/(m*l^2)) * (u(i) - b*x(2) - m*g*l*sin(x(1)));
    x(1) = x(1) + dt * dx(1);
    x(2) = x(2) + dt * dx(2);
    x_traj(:,i) = x;
end

end

function dX = pendulum_dynamics(t,X,u)
global m l g b
dX = zeros(2,1);
dX(1) = X(2);
dX(2) = (1/(m*l^2)) * (u - -b*X(2) - m*g*l*sin(X(1)));

end



function animation(theta_traj)


rodPivotPoint = [2 2]; %rectangular coordinates
rodLength = 1;
radius = .2;

theta0 = 0;
position = rodPivotPoint - (rodLength*[-sin(theta0) cos(theta0)]); %in rectangular coordinates

color_arr = [0,0,0];
% color_arr = 'red';
 
%Run simulation, all calculations are performed in cylindrical coordinates
    %Generate graphics, render pendulum
h = figure;
axesHandle = gca;
xlim(axesHandle, [(rodPivotPoint(1) - rodLength - radius) (rodPivotPoint(1) + rodLength + radius)] );
ylim(axesHandle, [(rodPivotPoint(2) - rodLength - radius) (rodPivotPoint(2) + rodLength + radius)] );

hold on
plot(rodPivotPoint(1),rodPivotPoint(2),'^'); %pendulum pivot
lineHandle = line([rodPivotPoint(1) position(1)],...
    [rodPivotPoint(2) position(2)],'linewidth',6,'color',color_arr); %pendulum rod
hold off
for i = 1:size(theta_traj,2)
    drawnow; %Forces MATLAB to render the pendulum
    position = rodPivotPoint - (rodLength*[-sin(theta_traj(i)) cos(theta_traj(i))]);
    set(lineHandle,'XData',[rodPivotPoint(1) position(1)],'YData',...
        [rodPivotPoint(2) position(2)],'linewidth',6,'color',color_arr);

    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 

    fn = strcat('iter', int2str(j), '.gif');
    % Write to the GIF File 
    if i == 1 
        imwrite(imind,cm,fn,'gif', 'Loopcount',inf); 
    else 
        imwrite(imind,cm,fn,'gif','WriteMode','append'); 
    end 
end

end

