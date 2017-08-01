function SolveWell_simple(welli,CurT)
global ActEle
global well Mat Fluid  
nFluid = Fluid.nfluid;
nProp = Fluid.nprop;
for i = 1 : well{welli}.nSch
    if CurT >= well{welli}.Sch(i).t0 && CurT <= well{welli}.Sch(i).t1
        Schindex = i;
    end
end
proptype =  well{welli}.Sch(Schindex).Prop;
fluidtype = well{welli}.Sch(Schindex).Fluid;
nPerf  = well{welli}.Sch(Schindex).nPf;
Perf = well{welli}.Sch(Schindex).Pf;
fdens = zeros(nFluid,1);
frefp = zeros(nFluid,1);
fcmp = zeros(nFluid,1);
propdens = zeros(nProp,1);

for i = 1 : nFluid
    name = Fluid.fluid{i}.name;
    fdens(i) = Fluid.fluid{i}.dens;
    frefp(i) = Fluid.fluid{i}.refp;
    fcmp(i) = Fluid.fluid{i}.compr;
    if strcmp(name,fluidtype) == 1
        fluidindex = i;
    end
end
%
for i = 1 : nProp
    name = Fluid.prop{i}.name;
    propdens(i) = Fluid.prop{i}.dens;
    if strcmp(name,proptype) == 1
        propindex = i;
    end
end

xfluid = zeros(nFluid,1);
xfluid(fluidindex) = 1;
cprop = zeros(nProp,1);
cprop(propindex) =  well{welli}.Sch(Schindex).PropFraction;
% m^3/s
rate = well{welli}.Sch(Schindex).ContrValue*0.159/60/nPerf;
[dens_sl0,~] = calcSLdens(Mat.Pp*1e6,nFluid,fdens,fcmp,frefp,xfluid,cprop,propdens);
mass_rate = dens_sl0*rate;
for i = 1 : nPerf
     
    numPf = Perf(i);
    numFrac = well{welli}.Perfindex(numPf,:);
    well{welli}.Sch(Schindex).PfQsl(i,1) = 0;
    well{welli}.Sch(Schindex).Pf_rate(i,1) =  0;
    well{welli}.Sch(Schindex).Pf_Q(i,1) =  0;
    well{welli}.Sch(Schindex).PfQsl(i,2) = 0;
    well{welli}.Sch(Schindex).Pf_rate(i,2) = 0;
    well{welli}.Sch(Schindex).Pf_Q(i,2) = 0;
   % if ActEle(numFrac(1)) > 0.1
        well{welli}.Sch(Schindex).PfQsl(i,1) = mass_rate/2;
      %  well{welli}.Sch(Schindex).Pf_rate(i,1) = veloc;
        well{welli}.Sch(Schindex).Pf_Q(i,1) = rate/2;
        fprintf('Transfer Mass = :%f\n',mass_rate);
  %  end
  %  if ActEle(numFrac(2)) > 0.1
        well{welli}.Sch(Schindex).PfQsl(i,2) = mass_rate/2;
      %  well{welli}.Sch(Schindex).Pf_rate(i,2) = veloc;
        well{welli}.Sch(Schindex).Pf_Q(i,2) = rate/2;
        fprintf('Transfer Mass = :%f\n',mass_rate);
   % end
    for mm = 1 : nProp
        well{welli}.Sch(Schindex).PfCp(i,mm) = cprop(mm);
    end
    for mm = 1 : nFluid
        well{welli}.Sch(Schindex).PfXf(i,mm) = xfluid(mm);
    end
end
end