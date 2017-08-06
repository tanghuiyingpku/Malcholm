function [Tau,Sn] = CalcPointStress_C_BC_local_3(num,DD)
global MaxEle;
global UseHeight;
 
global nAct CM_global GM_global IndexInv ;
%Stress Induced by Fracture deformation
CM1 = zeros(1,nAct*2);
CM2 = zeros(1,nAct*2);

for i = 1 : nAct
    CM1(i) = CM_global(IndexInv(num),IndexInv(i));
    CM1(i+nAct) = CM_global(IndexInv(num),IndexInv(i)+MaxEle);
    CM2(i) = CM_global(IndexInv(num)+MaxEle,IndexInv(i));
    CM2(i+nAct) = CM_global(IndexInv(num)+MaxEle,IndexInv(i)+MaxEle);
    if UseHeight > 0.1
        CM1(i) = CM1(i)*GM_global(IndexInv(num),IndexInv(i));
        CM1(i+nAct) = CM1(i+nAct)*GM_global(IndexInv(num),IndexInv(i));
        CM2(i) = CM2(i)*GM_global(IndexInv(num),IndexInv(i));
        CM2(i+nAct) = CM2(i+nAct)*GM_global(IndexInv(num),IndexInv(i));
    end
end

Tau = CM1*DD/1e6;
Sn =  CM2*DD/1e6;

clear CM1 CM2 DD;

end