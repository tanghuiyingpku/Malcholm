function [L1,w,p] = AnalyticalKGD(t)
Sh = 30.6897;
visco = 0.01;
rate = 10*0.159/60;
pr = 0.2;
E = 46.9101*1e9;
E = E/(1-pr^2);
h =66.0066;
L1 = 0.539*(rate^3*E/visco/h^3)^(1/6)*t.^(2/3);
w = 2.36*(visco*rate^3/E/h^3)^(1/6)*t.^(1/3);
p = Sh + 1.09*(E^2*visco)^(1/3)*t.^(-1/3)/1e6;
end