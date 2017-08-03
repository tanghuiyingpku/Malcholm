function [Kic,angle]= FindMaxTheta(a,b,KI1,KI2)
global Mat;
if abs(KI1) < 1e-3 && abs(KI2) < 1e-3
    Kic = 0;
    angle = 0;
    return;
end
a = -90;
b = -30;
[hasFind,c]= FindOpenTheta(a,b,KI1,KI2);
while hasFind < 0.1
    a = a+60;
    b = b+60;
    [hasFind,c]= FindOpenTheta(a,b,KI1,KI2);
end
if hasFind < 0.1
    %没有找到临界角度
    error('Fail to Find the critical Angle');
else
    k1 = 0.5*cosd(c/2)*(KI1*(1+cosd(c))-3*KI2*sind(c));
    k2 = 0.5*cosd(c/2)*(KI1*sind(c)+KI2*(3*cosd(c)-1));
    Kic = (k1/Mat.K1HF)^2+(k2/Mat.K2HF)^2;
    angle = c;
end
end