function [status, coord] = FindCrossPoint(La, Lb)
x1 = La(1);y1 = La(2);x2 = La(3);y2 = La(4);
x3 = Lb(1);y3 = Lb(2);x4 = Lb(3);y4 = Lb(4);
status = 0;
crossX = 0; 
crossY = 0;
coord = [crossX,crossY];
eps = 1.0E-8;
Marker = 0 ;
vector1 = [x2 - x1, y2 - y1];
vector2 = [x4 - x3, y4 - y3];
f = vector1(1) * vector2(2) - vector1(2) * vector2(1);
if abs(f) > eps
%     status  = 1 ;
    %K1 K2 гавтвх
    if abs(vector1(1)) < 1e-8
        Marker = 1;
        k2 = (y4-y3)/(x4-x3);
        b2 = y3 - k2*x3;
        crossX = x1;
        crossY = k2*crossX + b2;
    end
    if abs(vector2(1)) < 1e-8
        Marker = 1;
        k1 = (y2-y1)/(x2-x1);
        b1 = y1 - k1*x1;
        crossX = x3;
        crossY = k1*crossX + b1;
    end
    if Marker < eps
        k1 = (y2-y1)/(x2-x1);
        b1 = y1 - k1*x1;
        k2 = (y4-y3)/(x4-x3);
        b2 = y3 - k2*x3;
        crossX = (b1 - b2)/(k2-k1);
        crossY = k1*crossX + b1;
    end
else
    %ALMOST parallel
    status = 0;
    return;
end
coord = [crossX, crossY];
isInA = isDotIn(coord,La);
isInB = isDotIn(coord,Lb);
if isInA > 0.1 && isInB > 0.1
    status = 1;
end
end