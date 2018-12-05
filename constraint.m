function [c,ceq] = constraint(s, n, dt)

p = 4;
c = [];
% ceq = [s(1:4)-x0;
%        s((n-1)*4+1:n*4)-xf];
ceq = [];
for i = 1:n-1
  
    x = [s((i-1)*p+1); s((i-1)*p+2);s((i-1)*p+3);s((i-1)*p+4)];

    xn = [s(i*p+1); s(i*p+2);s(i*p+3);s(i*p+4)];

    u =  s(n*p+i);
    un = s(n*p+i+1);
    
    f = furataDynamics(x,u);
    fn = furataDynamics(xn,un);
    
    xhalf = (xn + x)/2 + dt*(f-fn)/8;
    uhalf = (u +  un)/2;
    fhalf = furataDynamics(xhalf,uhalf);
    ceq_ = x - xn + dt*(f+4*fhalf+fn)/6;
    ceq = [ceq;ceq_];
end

end