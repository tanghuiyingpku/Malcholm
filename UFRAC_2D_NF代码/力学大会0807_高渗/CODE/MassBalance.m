function detQ = MassBalance(CurT, e0,e,Schindex ,dens0,dens,h,dt,CL)
global AllEle_global  nAct Index IndexInv Tau_global
global well nwell
MassInj = 0;
for ii = 1 : nwell
    for jj = 1 :well{ii}.Sch(Schindex(ii)).nPf
        iele =  well{ii}.Sch(Schindex(ii)).Pf(jj);
        num = Index(well{ii}.Perfindex(iele,:));
        if min(num) > 1e-6
            MassInj = MassInj + well{ii}.Sch(Schindex(ii)).PfQsl(jj,1);
            MassInj = MassInj + well{ii}.Sch(Schindex(ii)).PfQsl(jj,2);
        end
    end
end
MassInj= MassInj*dt;
MassAccum = 0;
for i = 1 : nAct
    Li = AllEle_global(IndexInv(i),7);
    dmass = e(i)*h*Li*dens(i) - e0(i)*h*Li*dens0(i);
    MassAccum = MassAccum + dmass;
    MassAccum = MassAccum + 2*dt*CL*Li*h/sqrt(CurT +dt - Tau_global(IndexInv(i)))*dens(i);
end
detQ = (MassInj - MassAccum)/MassInj;
end