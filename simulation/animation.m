clear;
close all;
clc;

Lr = 1;
sita_degree = 30;
sita_rad = deg2rad(sita_degree);


p0 = [0;0;0];
d01 = [Lr;0;0];
p1 = rotz(sita_degree) * d01;



Lp = 2;
alpha_degree = 90;
d12 = [0;0;-Lp];
p2 = p1 + rotz(sita_degree) * rotx(alpha_degree) * d12;



% gather the coordinates and make them into x y z representaion
x = [p0(1), p1(1), p2(1)];
y = [p0(2), p1(2), p2(2)];
z = [p0(3), p1(3), p2(3)];


% plot the axis
plot3([0 1],[0 0],[0 0],'r','LineWidth',1.5);
hold on;
plot3([0 0],[0 1],[0 0],'g','LineWidth',1.5);
plot3([0 0],[0 0],[0 1],'b','LineWidth',1.5);
axis([-3,3,-3,3,-3,3]);


scatter3(x,y,z,'filled','LineWidth',5);
plot3(x,y,z,'black','LineWidth',2);
grid on;




function trajectory_animation(theta_traj)


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






