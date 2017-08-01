function [Gtotals,Gtotaln] = calcEnergy()
global   nAct ;
global   Mat DD_global  ;
global   AllEle_global  MaxEle PresF_global IndexInv ;
Pres = zeros(nAct,1);
DD = zeros(nAct*2,1);
InsituS = zeros(nAct*2,1);
AllEle = zeros(nAct,11);
% Project the in-situ stresses to fractures;
Sxx = Mat.Sxx*1e6;
Syy = Mat.Syy*1e6;
Sxy = Mat.Sxy*1e6;

for i = 1 : nAct
    Pres(i) = PresF_global(IndexInv(i))*1e6;
    DD(i) = DD_global(IndexInv(i));
    DD(i+nAct) = DD_global(IndexInv(i)+MaxEle);
    AllEle(i,:) = AllEle_global(IndexInv(i),:);
    %alpha = theta + 90;
    sinalp = AllEle(i,6);  %cosbet(i);
    cosalp = -AllEle(i,5); %-sinbet(i);
    InsituS(i) = (Sxx - Syy) * cosalp * sinalp - Sxy*(cosalp^2 - sinalp^2);
    InsituS(i + nAct) = Sxx * cosalp^2 + Syy * sinalp^2+2*Sxy*sinalp*cosalp;
end
Gtotals = 0;
Gtotaln = 0;

for i = 1 : nAct
    Gtotals = Gtotals + 0.5*AllEle(i,7)*(DD(i)*(0-InsituS(i)))*Mat.h;
    Gtotaln = Gtotaln + 0.5*AllEle(i,7)*(-Pres(i)-InsituS(i+nAct))*DD(i+nAct)*Mat.h;
end

end