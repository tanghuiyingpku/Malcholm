function resetValue()
global PresF_global Mat nAct stepL AllEle_global IndexInv;
PresF_global = PresF_global * 0 ;
L = stepL * 8 ;
rate = 10*0.159/60/2;
visco = 0.1;
for i = 1 : nAct
    PresF_global(IndexInv(i)) = 1.8*(16*visco*rate*(Mat.E*1e9)^3/pi/50^4*(L-abs(AllEle_global(IndexInv(i),9)-200)))^0.25-Mat.Sxx*1e6;
    PresF_global(IndexInv(i)) = PresF_global(IndexInv(i))/1e6;
end
CalcDD();
end