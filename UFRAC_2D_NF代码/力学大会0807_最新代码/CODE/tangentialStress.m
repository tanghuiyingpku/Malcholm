function tangentialStress()
theta = -90:1:90;
KI = 1;
KII = 2;
ts = 0.5*cosd(theta/2).*(KI*(1+cosd(theta)) - KII*3*sind(theta));
figure;
plot(theta,ts,'.')

end