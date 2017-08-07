function SolveWell_simple(welli,CurT)
global ActEle
global well Mat Fluid  DD_global PresF_global MaxEle
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
fvisco = zeros(Fluid.nfluid,1);

for i = 1 : nFluid
    fvisco(i) = Fluid.fluid{i}.visco;
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
rate = well{welli}.Sch(Schindex).ContrValue*0.159/60;
[dens_sl0,~] = calcSLdens(Mat.Pp*1e6,nFluid,fdens,fcmp,frefp,xfluid,cprop,propdens);
mass_rate = dens_sl0*rate;

kwf = zeros(nPerf,1);
Pf = zeros(nPerf,1);
for i = 1 : nPerf
    numPf = Perf(i);
    numFrac = well{welli}.Perfindex(numPf,:);
    well{welli}.Sch(Schindex).PfQsl(i,1) = 0;
    well{welli}.Sch(Schindex).Pf_rate(i,1) =  0;
    well{welli}.Sch(Schindex).Pf_Q(i,1) =  0;
    well{welli}.Sch(Schindex).PfQsl(i,2) = 0;
    well{welli}.Sch(Schindex).Pf_rate(i,2) = 0;
    well{welli}.Sch(Schindex).Pf_Q(i,2) = 0;
    if ActEle(numFrac(1)) > 0.1
        e = -DD_global(numFrac(1)+MaxEle);
        kwf(i) = e^3/12;
        Pf(i) = PresF_global(numFrac(1));
    end
    for mm = 1 : nProp
        well{welli}.Sch(Schindex).PfCp(i,mm) = cprop(mm);
    end
    for mm = 1 : nFluid
        well{welli}.Sch(Schindex).PfXf(i,mm) = xfluid(mm);
    end
end
% Solve for Pw (assume to be constant)
A = 0;
B = 0;
visco= dot(xfluid,fvisco);
Pw0 = -Mat.Sxx*1e6;
coeff = 1e-5;
if min(Pf) > 0.1
    [dens_sl,~] = calcSLdens(Pw0,nFluid,fdens,fcmp,frefp,xfluid,cprop,propdens);
    for i = 1 : nPerf
        A = A + 2*coeff*dens_sl/visco*kwf(i)*Mat.h;
        B = B + 2*coeff*dens_sl/visco*kwf(i)*Mat.h*Pf(i)*1e6;
    end
    Pw =( mass_rate + B)/A;
end
for i = 1 : nPerf
    numPf = Perf(i);
    numFrac = well{welli}.Perfindex(numPf,:);
    mass = 0;
    if ActEle(numFrac(1)) > 0.1
        mass = 1*coeff*dens_sl/visco*kwf(i)*Mat.h*(Pw - Pf(i)*1e6);
        well{welli}.Sch(Schindex).PfQsl(i,1) = mass;
      %  well{welli}.Sch(Schindex).Pf_rate(i,1) = veloc;
        well{welli}.Sch(Schindex).Pf_Q(i,1) = 1*coeff/visco*kwf(i)*Mat.h*(Pw - Pf(i)*1e6);
        fprintf('Transfer Mass = :%f\n',mass);
    end
    if ActEle(numFrac(2)) > 0.1
        well{welli}.Sch(Schindex).PfQsl(i,2) = 1*coeff*dens_sl/visco*kwf(i)*Mat.h*(Pw - Pf(i)*1e6);
      %  well{welli}.Sch(Schindex).Pf_rate(i,2) = veloc;
        well{welli}.Sch(Schindex).Pf_Q(i,2) = 1*coeff/visco*kwf(i)*Mat.h*(Pw - Pf(i)*1e6);
        fprintf('Transfer Mass = :%f\n',mass);
    end
    if mass < -0.1
        keyboard;
    end
end
% RHS  = 0 ;
% Jac = 0;
% Pw = Pw0;
% if sum(ActEle) > 0.1
%     %SI Unit
%     eps = 1e4;
%     while eps > 1e-4
%         for i = 1 : nPerf
%             [dens_sl,dens_dp] = calcSLdens(Pw,nFluid,fdens,fcmp,frefp,xfluid,cprop,propdens);
%             visco= dot(xfluid,fvisco);
%             RHS = RHS + 2*dens_sl/visco*kwf(i)*Mat.h*(Pw - Pf(i)*1e6);
%             Jac = Jac + 2*dens_sl/visco*kwf(i)*Mat.h*(1) + 2*dens_dp/visco*kwf(i)*Mat.h*(Pw - Pf(i)*1e6);
%         end
%         RHS = RHS - mass_rate;
%         dx = -Jac\RHS;
%         Pw = Pw + dx;
%         eps = norm(dx);
%     end
% end

    
    
    
end