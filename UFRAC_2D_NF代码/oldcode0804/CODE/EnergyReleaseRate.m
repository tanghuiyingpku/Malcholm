function EnergyReleaseRate()
% S = a11*KI + 2*a12*KI*K2 + a22*KII^2
close all;
theta = -90:1:90;
S = zeros(length(theta),1);
pr = 0.2;
k = 3 - 4*pr;
KI = 1;
KII = 2;
for i = 1 : length(theta)
    a11 = 1/16*(1+cosd(theta(i)))*(k - cosd(theta(i)));
    a12 = 1/16*sind(theta(i))*(2*cosd(theta(i)) - (k-1));
    a22 = 1/16*((k+1)*(1-cosd(theta(i))) + (1+cosd(theta(i)))*(3*cosd(theta(i))- 1));
    S(i) = a11*KI^2+2*a12*KI*KII+a22*KII^2;
end
figure
plot(theta,S,'.');