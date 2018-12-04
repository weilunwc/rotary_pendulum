clear;
close all;
clc;

Lr = 1;
sita_degree = 30;
sita_rad = deg2rad(sita_degree);


p0 = [0;0;0];
d01 = [Lr;0;0];
p1 = rotz(sita_degree) * d01



Lp = 2;
alpha_degree = 90;
d12 = [0;0;-Lp];
p2 = p1 + rotz(sita_degree) * rotx(alpha_degree) * d12



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



