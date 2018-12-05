clear all
close all
clc

syms th2 th1 thd2 thd1 thdd2 thdd1 V 'real'
%syms Dp Dr g Jp Jr km kt Lp Lr Mp Mr Rm 'real'
qube2_rotpen_param
x = [th1; th2; thd1; thd2];


l1 = Lr/2;
l2 = Lp/2; 
L1 = Lr;
L2 = Lp;
m1 = Mr;
m2 = Mp;
J1 = Jr;
J2 = Jp;
b1 = Dr;
b2 = Dp;

J1h = J1 + m1*l1^2;
J2h = J2 + m2*l2^2;
J0h = J1h+ m2*L1^2;

tau = km*(V-km*thd1)/Rm;

eom1 = thdd1*(J0h + J2h*sin(pi-th2)^2) + thdd2*(m2*L1*l2*cos(pi-th2)) - m2*L1*l2*sin(pi-th2)*thd2^2 + thd1*thd2*J2h*sin(2*(pi-th2))+b1*thd1 == tau;
eom2 = thdd1*m2*L1*l2*cos(pi-th2) + thdd2*J2h - 0.5*thd1^2*J2h*sin(2*(pi-th2))+ b2*thd2 + g*m2*l2*sin(pi-th2) == 0;

eom = [eom1, eom2];
sol = solve(eom, [thdd1 thdd2]);
thdd1_sol = sol.thdd1;
thdd2_sol = sol.thdd2;

f = [thd1, thd2, thdd1_sol, thdd2_sol];
A = jacobian(f, x);
B = jacobian(f, V);

matlabFunction([thd1;thd2;thdd1_sol;thdd2_sol],...
    'file','furataDynamics_z.m',...
    'vars',{[th1;th2;thd1;thd2],V},...
    'outputs',{'dyn'});

matlabFunction(A,B,...
    'file','linearizedfurata_z.m',...
    'vars',{[th1;th2;thd1;thd2],V},...
    'outputs',{'pfpx','pfpu'});
