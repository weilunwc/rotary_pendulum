function dS = Snomdynamics(t,s,Q,R,xnom,unom,time)
S2 = [s(1:4) s(5:8) s(9:12) s(13:16)];
S1 = [s(17); s(18); s(19); s(20)];
S0 = s(21);
[xn,un]= interped_xu(xnom,unom,t,time);
[A,B] = linearizedfurata_z(xn,0);
dS2 = -(Q-S2*B*inv(R)*B'*S2+S2*A+A'*S2);
dS1 = -(-2*Q*(xn)+(A'-S2*B*inv(R)*B')*S1+2*S2*B*(un));
dS0 = -(xn'*Q*xn-0.25*S1'*B*inv(R)*B'*S1+S1'*B*(un));
dS = [dS2(:);dS1(:);dS0];
end