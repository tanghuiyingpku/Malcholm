function Vs= CalcStoreV(nAllEle,AllEle,Pre,Dn,Mat)
%Calculate the stored volumn of fluid
Vs = 0;
for i = 1 : nAllEle
    if AllEle(i,10) < 1.1
        if Dn(i) < 0
            Dn(i) = 0;
        end
        b = Calc_Bw(Pre(i),1);
        %fprintf('Density at element %d is %f:\n',i,b);
        Vs = Vs + b*Dn(i)*AllEle(i,7)* Mat.h;
    end
end
end