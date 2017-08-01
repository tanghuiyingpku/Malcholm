function [KI1,KI2] = StressIntensF()
global nTip TipStates TipStatesInv;
global IndexInv nAct Mat;
global AllEle_global DD_global  MaxEle nActTip;
Dsi  = zeros(nAct,1);
Dni = Dsi;
for i = 1 : nAct
    Dsi(i) = DD_global(IndexInv(i));
    Dni(i) = DD_global(IndexInv(i) + MaxEle);
end
KI1 = zeros(nTip,1);
KI2 = zeros(nTip,1);
itip = 0;
for ii = 1: nAct
    if TipStates(IndexInv(ii)) > 0.1
        itip = itip +1;
        % if TipStates(IndexInv(ii)) < 998
        Ds= Dsi(ii);
        Dn= Dni(ii);
        if Dn > -1e-7
            Dn = 0;
          %  Ds = 0;
        end
        d = AllEle_global(IndexInv(ii),7)/2;
        if TipStates(IndexInv(ii)) < 998 &&  TipStates(IndexInv(ii)) > 0.01
            TipStates(IndexInv(ii))  = itip;
        end
        if TipStates(IndexInv(ii)) > 998
            TipStates(IndexInv(ii)) = itip + 999;
        end
        num = itip;
        TipStatesInv(itip) = IndexInv(ii);
        if  AllEle_global(IndexInv(ii),9) < 88
           % continue;
        end
        KI1(num) = -Mat.G*1e3/4/(1-Mat.miu)*(2*pi/d)^0.5*Dn;
        KI2(num) = -Mat.G*1e3/4/(1-Mat.miu)*(2*pi/d)^0.5*Ds;
    end
end
nActTip = itip;
%nTip  = itip;