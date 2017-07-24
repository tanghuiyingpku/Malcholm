function [Sh,SH] = PrinStress2(Sxx,Syy,Sxy)
[~,n] = size(Sxx);
Sh = zeros(n,n);
SH = Sh;
for i = 1 : n
    for j = 1 : n
        a = 1;
        b = -(Sxx(i,j) + Syy(i,j));
        c = Sxx(i,j) * Syy(i,j) - Sxy(i,j)^2;
        s1 = (-b+sqrt(b^2-4*a*c))/2/a;
        s2 = (-b-sqrt(b^2-4*a*c))/2/a;
        Sh(i,j) = min(s1,s2);
        SH(i,j) = max(s1,s2);
    end
end