function dS = SABdynamics(t,S,A,B,Q,R)
S = reshape(S,size(A));
ds = -(Q+A'*S+S*A-S*B*inv(R)*B'*S);
dS = ds(:);
end