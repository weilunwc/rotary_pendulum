%% RRT implementation
%% 
sita_range = [-4 4];
sita_dot_range = [-2 2];
start = [0;0];
goal = [-pi;0];
goal2 = [pi;0];
goal_bias = 0.1;

% Initialize tree
q_init.coord = [0;0];
q_init.parent = [];
q_init.id = 1;
q_init.control = [];
q_init.parent_id = 0;

iterations = 4000;
Tree.points = repmat(q_init,iterations+1,1);



for i = 1:iterations
    % goal bias
    if rand() <= goal_bias
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
     
    q.coord = new_coord;
    q.parent = nearest_coord;
    q.parent_id = nearest_id;
    q.id = i+1;
    q.control = best_u;
    Tree.points(i+1) = q;
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
