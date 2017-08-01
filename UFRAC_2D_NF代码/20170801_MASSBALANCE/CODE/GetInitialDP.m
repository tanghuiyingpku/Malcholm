function [DD,P] = GetInitialDP()
global   nAct ;
global   DD_global  ;
global   MaxEle PresF_global IndexInv ;
P = zeros(nAct,1);
DD = zeros(nAct*2,1);

for i = 1 : nAct
    P(i) = PresF_global(IndexInv(i));
    DD(i) = DD_global(IndexInv(i));
    DD(i+nAct) = DD_global(IndexInv(i)+MaxEle);
end
end