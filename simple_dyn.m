function dx = simple_dyn(x,u)
A = [0 0 1 0;
     0 0 0 1;
     0 149.2751 -0.0104 0;
     0 261.6091 -0.0103 0];
B = [0;0;49.7275; 49.1493];
dx = A*x+B*u;
end