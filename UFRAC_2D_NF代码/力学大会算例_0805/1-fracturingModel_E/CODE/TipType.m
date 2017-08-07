function [type,vec] = TipType(ne,AllEle,ConnList)
global Tipcoordinate;
p0 = [AllEle(ne,8),AllEle(ne,9)];
p1 = Tipcoordinate(ne,:);
dis = p0-p1;
vec = dis/sqrt(dis(1)^2+dis(2)^2);
vec0 = [AllEle(ne,6) AllEle(ne,5)];
if dot(vec,vec0) < -1e-8
    %ÄæÏò
    type = 1;
else
    type = 2;
end

end