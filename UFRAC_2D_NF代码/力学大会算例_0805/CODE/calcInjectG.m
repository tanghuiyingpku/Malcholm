function  calcInjectG(CurT,dt)
global well  nwell G_inj Index Mat
Sxx = Mat.Sxx;
Schindex = zeros(nwell,1); 
for welli = 1 : nwell
    for i = 1 : well{welli}.nSch
        if CurT >= well{welli}.Sch(i).t0 && CurT <= well{welli}.Sch(i).t1
            Schindex(welli) = i;
        end
    end
end
for ii = 1 : nwell
    for jj = 1 :well{ii}.Sch(Schindex(ii)).nPf
        iele =  well{ii}.Sch(Schindex(ii)).Pf(jj);
        num = Index(well{ii}.Perfindex(iele,:));
        if min(num) > 1e-6 && well{ii}.Sch(Schindex(ii)).Pf_Q(jj,1) > 1e-11
            G_inj = G_inj + well{ii}.Sch(Schindex(ii)).Pf_Q(jj,1) * 2 * (well{ii}.Sch(Schindex(ii)).Pres(jj)+Sxx) * dt* 1e6;
        end
    end
end
end
