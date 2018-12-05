clear all
close all
clc

syms theta theta_dot theta_ddot alpha alpha_dot alpha_ddot tau 'real'
%syms Dp Dr g Jp Jr km kt Lp Lr Mp Mr Rm 'real'
qube2_rotpen_param
s =[theta;alpha];
ds = [theta_dot; alpha_dot];
dds = [theta_ddot; alpha_ddot];

% i = sym([1;0;0]);
% j = sym([0;1;0]);
% k = sym([0;0;1]);
% 
% vc1 = -j*Lr/2*theta_dot;
% w2 = i
% rr = Lr*(cos(theta)*i+sin(theta)*j);
% rc = rr/2;
% rp = rr-Lp/2*sin(alpha)*cos(theta)*i-Lp/2*sin(alpha)*sin(theta)*j-Lp/2*cos(alpha)*k;
% 
% s =[theta;alpha];
% ds = [theta_dot; alpha_dot];
% dds = [theta_ddot; alpha_ddot];
% vr = jacobian(rc,s)*ds;
% vp = jacobian(rp,s)*ds;
% %Trv = Mr*dot(vr,vr)/2;
% Trw = Jr*theta_dot^2/2;
% Tpv = Mp*dot(vp,vp)/2;
% Tpw = Jp*alpha_dot^2/2;
% 
% T = Trw+Tpv+Tpw;
% U = -Mp*g*Lp/2*cos(alpha);
% L = T-U;
% 
% JL = jacobian(L,ds);
% 
% eom = jacobian(JL,[s;ds])*[ds;dds] == (jacobian(L,s)' + [tau-Dr*theta_dot;-Dp*alpha_dot]);
% sol = solve(eom,dds);
% 
% tddot = sol.theta_ddot;
% addot = sol.alpha_ddot

% eom1 = (Mr*Lr^2+0.25*Mp*Lp^2-0.25*Mp*Lp^2*cos(alpha)^2+Jr)*theta_ddot...
%       -(0.5*Mp*Lp*Lr*cos(alpha))*alpha_ddot+(0.5*Mp*Lp^2*sin(alpha)*cos(alpha))*theta_dot*alpha_dot...
%       +(0.5*Mp*Lp*Lr*sin(alpha))*alpha_dot^2 == tau - Dr*theta_dot;
% eom2 = 0.5*Mp*Lp*Lr*cos(alpha)*theta_ddot+(Jp+0.25*Mp*Lp^2)*alpha_ddot-0.25*Mp*Lp^2*cos(alpha)*sin(alpha)*theta_dot^2 ...
%       +0.5*Mp*Lp*g*sin(alpha) == -Dp*alpha_dot;

eom1 = (Mp*Lr^2+0.25*Mr*Lr^2-0.25*Mp*Lp^2*cos(alpha)^2+Jr)*theta_ddot...
       -(0.5*Mp*Lp*Lr*cos(alpha))*alpha_ddot+(0.5*Mp*Lp^2*sin(alpha)*cos(alpha))*theta_dot*alpha_dot...
       +(0.5*Mp*Lp*Lr*sin(alpha))*alpha_dot^2 == tau - Dr*theta_dot;
eom2 = 0.5*Mp*Lp*Lr*cos(alpha)*theta_ddot+(Jp+0.25*Mp*Lp^2)*alpha_ddot-0.25*Mp*Lp^2*cos(alpha)*sin(alpha)*theta_dot^2 ...
       +0.5*Mp*Lp*g*sin(alpha) == -Dp*alpha_dot;

eom1 = (Jr+0.25*Mr*Lr^2+Mp*Lr^2)*theta_ddot...
      +
eom2 = 0.5*Mp*Lp*Lr*cos(alpha)*theta_ddot+(Jp+0.25*Mp*Lp^2)*alpha_ddot-0.25*Mp*Lp^2*cos(alpha)*sin(alpha)*theta_dot^2 ...
      +0.5*Mp*Lp*g*sin(alpha) == -Dp*alpha_dot;

%eom1 = (Mp*Lr^2+Jr)*theta_ddot-Mp*Lp*Lr*alpha_ddot/2 == tau-Dr*theta_dot;
%eom2 = Mp*Lp*Lr*theta_ddot/2 + (Jp+Mp*Lp^2/4)*alpha_ddot+Mp*Lp*g*alpha/2 == -Dp*alpha_dot;
eom = [eom1, eom2];
sol = solve(eom, [theta_ddot alpha_ddot]);
tddot = sol.theta_ddot;
addot = sol.alpha_ddot;
A3 = jacobian(tddot,[s;ds])
A4 = jacobian(addot,[s;ds])
B3 = jacobian(tddot,tau)
B4 = jacobian(addot,tau)

A = [0 0 1 0;
     0 0 0 1;
     A3;
     A4];
B = [0 ;
     0 ;
     B3;
     B4];
%Jt  = Jp*Mp*Lr^2+Jr*Jp+Jr*Mp*Lp^2/4;
%subs(tddot, [theta,alpha,theta_dot,alpha_dot],[0,0,0,0])
%tddot = (-(Jp+Mp*Lp^2/4 
% f = [theta_dot; alpha_dot;tddot;addot];
% pfpx = jacobian(f,[s;ds]);   
% pfpu = jacobian(f,tau);

%A = subs(pfpx, [theta, alpha, theta_dot, alpha_dot] , [theta_query, alpha_query, theta_dot_query, alpha_dot_query]);


matlabFunction([theta_dot;alpha_dot;tddot;addot],...
    'file','furataDynamics.m',...
    'vars',{[theta; alpha; theta_dot; alpha_dot],tau},...
    'outputs',{'dyn'});

matlabFunction(A,B,...
    'file','linearizedfurata.m',...
    'vars',{[theta; alpha; theta_dot; alpha_dot],tau},...
    'outputs',{'pfpx','pfpu'});