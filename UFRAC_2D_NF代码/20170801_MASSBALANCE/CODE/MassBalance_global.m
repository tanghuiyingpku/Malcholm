function  err = MassBalance_global(CurT)
global AllEle_global nAct e_global einit_global IndexInv
massStore = 0;
massInj = 10*0.159*CurT/60;

for i  = 1 : nAct
    massStore = massStore + e_global(IndexInv(i))*50*AllEle_global(IndexInv(i),7);
    massStore = massStore - einit_global(IndexInv(i))*50*AllEle_global(IndexInv(i),7);
end

err = (massInj - massStore)/massInj;
end