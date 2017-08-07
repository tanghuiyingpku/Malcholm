function [Kic,angle]= FindMaxTheta(a,b,KI1,KI2)
global Mat;
Kic = 0;
angle = 0;
temp = a : 0.1: b;
Kni = zeros(length(temp),1);
num = 0;
for theta0 = a : 0.1: b
    num = num + 1;
    k1 = 0.5*cosd(theta0/2)*(KI1*(1+cosd(theta0))-3*KI2*sind(theta0));
    k2 = 0.5*cosd(theta0/2)*(KI1*sind(theta0)+KI2*(3*cosd(theta0)-1));
    Ki = (k1/Mat.K1HF)^2+(k2/Mat.K2HF)^2;
    Kni(num) = Ki;
    if Ki > Kic && k1 > 1e-10
        Kic = Ki;
        angle = theta0;
    end
end
end