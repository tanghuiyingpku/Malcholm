function L= CalculateDis(Point1,Point2)
[nPoint,~] = size(Point1);
L = zeros(nPoint,1);
for i = 1 : nPoint
    L(i) = sqrt((Point2(i,1)-Point1(i,1))^2+(Point2(i,2)-Point1(i,2))^2);
end
end