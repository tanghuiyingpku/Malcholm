function detQ = FlowBalance(nAllEle,AllEle,DD0,DD,Pres,ConnList,dt)
detQ = zeros(nAllEle,1);
Pres = Pres*1e6;
for i = 1 : nAllEle
    numL = ConnList(i,3);
    numR = ConnList(i,4);
    ei = DD(i);
    ki = ei^2/12;
    Li = AllEle(i,7);
    if numL > 0.1
        ej = DD(numL);
        kj = ej^2/12;
        revmiu = Calc_miu(Pres(i),1);
        Lj = AllEle(numL,7);
        kij = ki*kj*(Li+Lj)/(ki*Li+kj*Lj);
        A = (DD(i)+DD(numL))/2;
        dx = 0.5*(AllEle(i,7)+AllEle(numL,7));
        detQ(i) = detQ(i) + revmiu*kij*A*(Pres(numL)-Pres(i))/dx;
    end
    if numR > 0.1
        ej = DD(numR);
        kj = ej^2/12;
        revmiu = Calc_miu(Pres(i),1);
        Lj = AllEle(numR,7);
        kij = ki*kj*(Li+Lj)/(ki*Li+kj*Lj);
        A = (DD(i)+DD(numR))/2;
        dx = 0.5*(AllEle(i,7)+AllEle(numR,7));
        detQ(i) = detQ(i) + revmiu*kij*A*(Pres(numR)-Pres(i))/dx;
    end
    detQ(i) = detQ(i) - (DD(i) - DD0(i))*AllEle(i,7)/dt;
end
