clear all
close all
clc

syms z1 z2 z3 z4 'real'
zd = [0;0;0;0];
[A,B] = linearizedfurata_z(zd,0);
Q = 10*eye(4);
R = 1;
[K,S] = lqr(A,B,Q,R);
z = [z1;z2;z3;z4];
f = furataDynamics_z(z,-K*z);
f_taylor = taylor(f,z,'ExpansionPoint',zd,'Order', 4);
J = z'*S*z;
Jdot = 2*z'*S*f_taylor;

rho = 1
INFO.feasratio = 1
while INFO.feasratio>0
rho = rho+1;
prog = sosprogram(z);
[prog,h] = sossosvar(prog,z);
expr =-(Jdot+h*(rho - J))-0.0001*sum(z.*z);
prog = sosineq(prog,expr);
solver_opt.solver = 'sedumi';
[prog,INFO] = sossolve(prog,solver_opt);
h_solved = sosgetsol(prog,h);
J_solved = sosgetsol(prog,J);
double(subs(J_solved, z,zd))
rho
end

