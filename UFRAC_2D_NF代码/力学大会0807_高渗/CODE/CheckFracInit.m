function isInit = CheckFracInit(Pw,iwell,CurT)
global Mat well  ActEle
isInit = 0;
nG = well{1}.nGrid;
Sh = min(abs(Mat.Sxx),abs(Mat.Syy));
SH = abs(Mat.Syy);
% alpha = 1;
CritP = SH;%(3*Sh - SH + Mat.tens - alpha*Mat.Pp*(1-2*Mat.miu)/(1-Mat.miu))/(1+Mat.poro-alpha*(1-2*Mat.miu)/(1-Mat.miu));
nSch = well{iwell}.nSch;
for j = 1 : nSch
    if CurT >= well{iwell}.Sch(j).t0 && CurT <= well{iwell}.Sch(j).t1
        for i = 1 : well{iwell}.Sch(j).nPf
            rate = well{iwell}.Sch(j).ContrValue*0.159/60/well{iwell}.Sch(j).nPf;
            Perfi = well{iwell}.Sch(j).Pf(i);
            Pwf = Pw;
            well{iwell}.Sch(j).Pf_Q(i,1) = rate/2;
            well{iwell}.Sch(j).Pf_Q(i,2) = rate/2;
            [dens_sl0,CP]= CalcWellDens(iwell,j);
            well{iwell}.Sch(j).PfCp(i,:) = CP;
            well{iwell}.Sch(j).PfQsl(i,1) = rate*dens_sl0/2;
            well{iwell}.Sch(j).PfQsl(i,2) = rate*dens_sl0/2;
            if Pwf > CritP
                Elei = well{iwell}.Perfindex(Perfi);
                if ActEle(Elei) < 0.1
                    disp('*****************************');
                    fprintf('Fracture Initiation\n');
                    disp(' ');
                    disp('*****************************');
                    Precondition(CurT,Elei,rate);
                    isInit = 1;
                end
            end
        end
    end
end

end