function dyn = furataDynamics_z(in1,V)
%FURATADYNAMICS_Z
%    DYN = FURATADYNAMICS_Z(IN1,V)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    04-Dec-2018 20:10:39

th2 = in1(2,:);
thd1 = in1(3,:);
thd2 = in1(4,:);
t2 = cos(th2);
t3 = sin(th2);
t4 = th2.*2.0;
t5 = sin(t4);
t6 = t2.^2;
t7 = t3.^2;
t8 = t7.*9.42319734012454e36;
t9 = t6.*-9.20532730968088e36+t8+2.84683270503094e37;
t10 = 1.0./t9;
t11 = thd1.^2;
t12 = thd2.^2;
dyn = [thd1;thd2;t10.*(V.*-2.831319433965669e36+thd1.*9.683112464162587e35+t2.*t3.*8.499224556044183e36-t3.*t12.*7.450900222424055e34+t2.*thd2.*2.798397114966068e35+t2.*t5.*t11.*3.725450111212027e34-t5.*thd1.*thd2.*7.538557872099632e34).*-1.25e2;t10.*(t3.*2.029609653804761e43+thd2.*6.682555287559105e41-V.*t2.*2.18624774606724e42+t3.*t7.*6.718137057167277e42+t5.*t11.*8.896352203221686e40+t2.*thd1.*7.476967291549962e41+t7.*thd2.*2.211968307785679e41-t2.*t3.*t12.*5.75332956855055e40+t5.*t7.*t11.*2.944749168788919e40-t2.*t5.*thd1.*thd2.*5.821015798768792e40).*(-1.6e-4)];