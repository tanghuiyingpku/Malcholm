function Vfl = CalcFluidV(nEle,dens_slv,Vel_sl,Vprop_x,Cp,propdens)
Vfl = zeros(nEle,4);
nprop = length(propdens);
for i = 1 : nEle
    Vfl(i,:) = dens_slv(i)*Vel_sl(i,:);
    for j = 1 : nprop
        Vfl(i,:) = Vfl(i,:) - Cp(j)*propdens(j)*Vprop_x(j,:);
    end
    temp = dot(Cp(i,:),propdens);
    Vfl(i,:) = Vfl(i,:)/(dens_slv(i) - temp);
end
end