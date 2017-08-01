function [dens_sl0,cprop]= CalcWellDens(welli,Schindex)
global well Mat Fluid
nFluid = Fluid.nfluid;
nProp = Fluid.nprop;
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
[dens_sl0,~] = calcSLdens(Mat.Pp*1e6,nFluid,fdens,fcmp,frefp,xfluid,cprop,propdens);
end