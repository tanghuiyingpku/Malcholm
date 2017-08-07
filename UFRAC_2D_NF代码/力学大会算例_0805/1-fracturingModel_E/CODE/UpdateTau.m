function UpdateTau(nActOld, nAct,CurT)
global Tau_global;
global IndexInv;
for i = nActOld+1:nAct
    Tau_global(IndexInv(i)) = CurT;
end
end