function [KI1,KI2] = StressIntensF_local(DD)
global nTip TipStates TipStatesInv;
global IndexInv nAct Mat;
global AllEle_global   nActTip;
Dsi  =DD(1:nAct);
Dni = DD(nAct+1:2*nAct);
KI1 = zeros(nTip,1);
KI2 = zeros(nTip,1);
itip = 0;
for ii = 1: nAct
    if TipStates(IndexInv(ii)) > 0.1
        itip = itip +1;
        % if TipStates(IndexInv(ii)) < 998
        if Dni(ii) > -1e-7
            Dni(ii) =0;
        end
        Ds= Dsi(ii);
        Dn= Dni(ii);
        d = AllEle_global(IndexInv(ii),7)/2;
        if TipStates(IndexInv(ii)) < 998 &&  TipStates(IndexInv(ii)) > 0.01
            TipStates(IndexInv(ii))  = itip;
        end
        num = itip;
        TipStatesInv(itip) = IndexInv(ii);
        KI1(num) = -Mat.G*1e3/4/(1-Mat.miu)*(2*pi/d)^0.5*Dn;
        KI2(num) = -Mat.G*1e3/4/(1-Mat.miu)*(2*pi/d)^0.5*Ds;
    end
end
nActTip = itip;