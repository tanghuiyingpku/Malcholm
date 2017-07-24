function Frac2 = RefineFrac(nEle,Frac,refineN)
Frac2 = zeros(nEle*refineN,4);
for i = 1 : nEle
    xx  = linspace(Frac(i,1)*1.1,Frac(i,3)*0.9,refineN+1);
    yy  = linspace(Frac(i,2)*1.1,Frac(i,4)*0.9,refineN+1);
    for j = 1 : refineN
        Frac2((i-1)*refineN+j,2) = yy(j);
        Frac2((i-1)*refineN+j,4) = yy(j+1);
        Frac2((i-1)*refineN+j,1) = xx(j);
        Frac2((i-1)*refineN+j,3) = xx(j+1);
    end
end