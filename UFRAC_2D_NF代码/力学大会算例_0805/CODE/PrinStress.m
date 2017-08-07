function [Sh,SH] = PrinStress(Sxx,Syy,Sxy)
n = length(Sxx);
Sh = zeros(n,1);
SH = Sh;
for i = 1 : n
    a = 1;
    b = -(Sxx(i) + Syy(i));
    c = Sxx(i) * Syy(i) - Sxy(i)^2;
    s1 = (-b+sqrt(b^2-4*a*c))/2/a;
    s2 = (-b-sqrt(b^2-4*a*c))/2/a;
    Sh(i) = min(s1,s2);
    SH(i) = max(s1,s2);
end
end