function Jac2 = DumpJacob(Jac)
 Jac2 = Jac;
[nx,ny] = size(Jac);
for i = 1 : nx
    for j = 1 : ny
        if abs(Jac(i,j))/abs(Jac(i,i))< 1e-6
            Jac2(i,j) = 0;
        end
    end
end
end