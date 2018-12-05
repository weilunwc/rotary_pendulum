function [A,B] = linearized_furata_top()

qube2_rotpen_param

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

A31 = 0;
A32 = g*m2^2*l2^2*L1/(J0h*J2h-m2^2*L1^2*l2^2);
A33 = -b1*J2h/(J0h*J2h-m2^2*L1^2*l2^2);
A34 = -b2*m2*l2*L1/(J0h*J2h-m2^2*L1^2*l2^2);
A41 = 0;
A42 = g*m2*l2*J0h/(J0h*J2h-m2^2*L1^2*l2^2);
A43 = -b1*m2*l2*L1/(J0h*J2h-m2^2*L1^2*l2^2);
A44 = -b2*J0h/(J0h*J2h-m2^2*L1^2*l2^2);
B31 = J2h/(J0h*J2h-m2^2*L1^2*l2^2);
B41 = m2*L1*l2/(J0h*J2h-m2^2*L1^2*l2^2);

A = [0 0 1 0 ;
     0 0 0 1;
     A31 A32 A33 A34;
     A41 A42 A43 A44];
B = [0; 0; B31; B41];