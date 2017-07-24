function [Pw,Xf,Cp] = SolveWell(Pw0,Xf0,Cp0,welli,dt,CurT)
global PresF_global  DD_global  ActEle
global well Mat Fluid MaxEle stepL
nG = well{welli}.nGrid;
dx = well{welli}.dx;
nFluid = Fluid.nfluid;
nProp = Fluid.nprop;
fvisco = zeros(nFluid,1);
fdens = zeros(nFluid,1);
frefp = zeros(nFluid,1);
fcmp = zeros(nFluid,1);
propdens = zeros(nProp,1);
charaP = 1e6;%Characteristic Pressure
%Picard Iteration Speed
alphaR = 0;%0.1;

wellD = 1e-1;
Kw = 1e-3*pi/32*wellD^2;
Aw = pi/4*wellD^2;
%
for i = 1 : well{welli}.nSch
    if CurT >= well{welli}.Sch(i).t0 && CurT <= well{welli}.Sch(i).t1
        Schindex = i;
    end
end
proptype =  well{welli}.Sch(Schindex).Prop;
fluidtype = well{welli}.Sch(Schindex).Fluid;
nPerf  = well{welli}.Sch(Schindex).nPf;
Perf = well{welli}.Sch(Schindex).Pf;
% For viscosity calculation
cmax = 0.6;
n = 1.3;
for i = 1 : nFluid
    name = Fluid.fluid{i}.name;
    fvisco(i) = Fluid.fluid{i}.visco;
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
% m/s
rate = well{welli}.Sch(Schindex).ContrValue*0.159/60/Aw;
[dens_sl0,~] = calcSLdens(Mat.Pp*1e6,nFluid,fdens,fcmp,frefp,xfluid,cprop,propdens);
dens_flinj = calcSLdens(Mat.Pp*1e6,nFluid,fdens,fcmp,frefp,xfluid,cprop*0,propdens);

mass_rate = dens_sl0*rate;
% First Time Step
if CurT < 1e-6
    Xf0((fluidindex-1)*nG+1:fluidindex*nG) = 1;
    Cp0((propindex-1)*nG+1:propindex*nG) = cprop(propindex);
end
Xf = Xf0;
Cp = Cp0;

% Main Loop
eps2 = 1e6;
tol = 1e-4;
Pw = Pw0;
RHS0 = zeros(nG,1);
Vel = zeros(nG,1);
Jac = zeros(nG,nG);
MatCp = zeros(nG*nProp,nG*nProp);
RhsCp = zeros(nG*nProp,1);
MatXf = zeros(nG*nFluid,nG*nFluid);
RhsXf = zeros(nG*nFluid,1);
iIter2 = 0;
visco_slv = zeros(nG,1);
dens_slv = zeros(nG,1);
dens_slv0 = zeros(nG,1);
dens_dpv = zeros(nG,1);
densfv = zeros(nG,nFluid);
densf0v = zeros(nG,nFluid);
while eps2 > tol
    eps = 1e6;
    iIter2 = iIter2+ 1;
    if iIter2 > 30
        iIter2 = 1;
        dt = dt /2;
        Xf = Xf0;
        Cp = Cp0;
        Pw = Pw0;
    end
    iIter1 = 0;
    while eps > tol
        iIter1 = iIter1 + 1;
        RHS0 = RHS0*0;
        Jac = Jac*0;
        Vel = Vel*0;
        for i = 1 : nG
            xfluid = Xf(i:nG:(nFluid-1)*nG+i);
            xfluid0 = Xf0(i:nG:(nFluid-1)*nG+i);
            cprop =  Cp(i:nG:(nProp-1)*nG+i);
            cprop0 =  Cp0(i:nG:(nProp-1)*nG+i);
            cprops = sum(cprop);
            visco_flv= dot(xfluid,fvisco);
            visco_slv(i) = visco_flv*(1-cprops/cmax)^-n;
            for mm = 1 : nFluid
                densf0v(i,mm) = CalcFLDens(mm,Pw0(i));
                densfv(i,mm) = CalcFLDens(mm,Pw(i));
            end
            % Original RHS
            [dens_slv(i),dens_dpv(i)] = calcSLdens(Pw(i)*charaP,nFluid,fdens,fcmp,frefp,xfluid,cprop,propdens);
            [dens_slv0(i),~] = calcSLdens(Pw0(i)*charaP,nFluid,fdens,fcmp,frefp,xfluid0,cprop0,propdens);
            dens_dpv(i) = dens_dpv(i) * charaP;
        end
        % Flux from well to fracture
        for i = 1 : nPerf
            numPf = Perf(i);
            numFrac = well{welli}.Perfindex(numPf,:);
            numWell = well{welli}.PerfWellindex(numPf);
            well{welli}.Sch(Schindex).PfQsl(i,1) = 0;
            well{welli}.Sch(Schindex).Pf_rate(i,1) =  0;
            well{welli}.Sch(Schindex).Pf_Q(i,1) =  0;
            well{welli}.Sch(Schindex).PfQsl(i,2) = 0;
            well{welli}.Sch(Schindex).Pf_rate(i,2) = 0;
            well{welli}.Sch(Schindex).Pf_Q(i,2) = 0;
            if ActEle(numFrac(1)) > 0.1
                Pfrac = PresF_global(numFrac(1));
                w = -DD_global(MaxEle+numFrac(1));
                Kwf = w^3/12*Mat.h;
                coeff = 1e-2;
                veloc =coeff * w^2/12/visco_slv(numWell)*(Pw(numWell) - Pfrac)*1e6/stepL;
                massWF =  coeff * dens_slv(numWell)*Kwf/visco_slv(numWell)*(Pw(numWell) - Pfrac)*1e6/stepL;
                RHS0(numWell) = RHS0(numWell) - massWF;
                well{welli}.Sch(Schindex).PfQsl(i,1) = massWF;
                well{welli}.Sch(Schindex).Pf_rate(i,1) = veloc;
                well{welli}.Sch(Schindex).Pf_Q(i,1) = veloc * Mat.h*w;
                fprintf('Transfer Mass = :%f\n',massWF);
            end
            if ActEle(numFrac(2)) > 0.1
                Pfrac = PresF_global(numFrac(2));
                w = -DD_global(MaxEle+numFrac(2));
                Kwf = w^3/12*Mat.h;
                coeff = 1e-2;
                veloc =coeff * w^2/12/visco_slv(numWell)*(Pw(numWell) - Pfrac)*1e6/stepL;
                massWF =  coeff * dens_slv(numWell)*Kwf/visco_slv(numWell)*(Pw(numWell) - Pfrac)*1e6/stepL;
                RHS0(numWell) = RHS0(numWell) - massWF;
                well{welli}.Sch(Schindex).PfQsl(i,2) = massWF;
                well{welli}.Sch(Schindex).Pf_rate(i,2) = veloc;
                well{welli}.Sch(Schindex).Pf_Q(i,2) = veloc * Mat.h*w;
                fprintf('Transfer Mass = :%f\n',massWF);
            end
        end
        for i = 1 : nG
            %Injection Rate
            % Flow Term
            if i == 1
                RHS0(i) = RHS0(i) + mass_rate*Aw;
                Vel(i) = rate;
                if Pw(i+1) > Pw(i)
                    mm = i+ 1;
                    visco_sl = visco_slv(mm);
                    dens_slj = dens_slv(mm);
                    dens_dpj = dens_dpv(mm);
                    RHS0(i) = RHS0(i) + dens_slj*Kw/visco_sl*Aw*(Pw(i+1) - Pw(i))/dx*charaP;
                    Jac(i,i+1) = Jac(i,i+1) +  dens_slj*Kw/visco_sl*Aw/dx*charaP;
                    Jac(i,i+1) = Jac(i,i+1) +  dens_dpj*Kw/visco_sl*Aw*(Pw(i+1) - Pw(i))/dx*charaP;
                    Jac(i,i)  =  Jac(i,i) - dens_slj*Kw/visco_sl*Aw/dx*charaP;
                else
                    mm = i;
                    visco_sl = visco_slv(mm);
                    dens_sl = dens_slv(mm);
                    dens_dp = dens_dpv(mm);
                    RHS0(i) = RHS0(i) + dens_sl*Kw/visco_sl*Aw*(Pw(i+1) - Pw(i))/dx*charaP;
                    Jac(i,i+1) = Jac(i,i+1) +  dens_sl*Kw/visco_sl*Aw/dx*charaP;
                    Jac(i,i)  =  Jac(i,i) - dens_sl*Kw/visco_sl*Aw/dx*charaP + dens_dp*Kw/visco_sl*Aw*(Pw(i+1) - Pw(i))/dx*charaP;
                end
            end
            if i == nG
                if Pw(i-1) > Pw(i)
                    mm =  i- 1;
                    visco_sl = visco_slv(mm);
                    dens_slj = dens_slv(mm);
                    dens_dpj = dens_dpv(mm);
                    RHS0(i) = RHS0(i) + dens_slj*Kw/visco_sl*Aw*(Pw(i-1) - Pw(i))/dx*charaP;
                    Jac(i,i-1) = Jac(i,i-1) +  dens_slj*Kw/visco_sl*Aw/dx*charaP;
                    Jac(i,i-1) = Jac(i,i-1) +  dens_dpj*Kw/visco_sl*Aw*(Pw(i-1) - Pw(i))/dx*charaP;
                    Jac(i,i)  =  Jac(i,i) - dens_slj*Kw/visco_sl*Aw/dx*charaP;
                else
                    mm = i;
                    visco_sl = visco_slv(mm);
                    dens_sl = dens_slv(mm);
                    dens_dp = dens_dpv(mm);
                    RHS0(i) = RHS0(i) + dens_sl*Kw/visco_sl*Aw*(Pw(i-1) - Pw(i))/dx*charaP;
                    Jac(i,i-1) = Jac(i,i-1) +  dens_sl*Kw/visco_sl*Aw/dx*charaP;
                    Jac(i,i)  =  Jac(i,i) - dens_sl*Kw/visco_sl*Aw/dx*charaP + dens_dp*Kw/visco_sl*Aw*(Pw(i-1) - Pw(i))/dx*charaP;
                end
                Vel(i) = Kw/visco_sl*(Pw(i-1) - Pw(i))/dx*charaP;
            end
            if i > 1.1 && i < nG-0.1
                
                if Pw(i+1) > Pw(i)
                    mm = i+ 1;
                    visco_sl = visco_slv(mm);
                    dens_slj = dens_slv(mm);
                    dens_dpj = dens_dpv(mm);
                    RHS0(i) = RHS0(i) + dens_slj*Kw/visco_sl*Aw*(Pw(i+1) - Pw(i))/dx*charaP;
                    Jac(i,i+1) = Jac(i,i+1) +  dens_slj*Kw/visco_sl*Aw/dx*charaP;
                    Jac(i,i+1) = Jac(i,i+1) +  dens_dpj*Kw/visco_sl*Aw*(Pw(i+1) - Pw(i))/dx*charaP;
                    Jac(i,i)  =  Jac(i,i) - dens_slj*Kw/visco_sl*Aw/dx*charaP;
                else
                    mm = i;
                    visco_sl = visco_slv(mm);
                    dens_sl = dens_slv(mm);
                    dens_dp = dens_dpv(mm);
                    RHS0(i) = RHS0(i) + dens_sl*Kw/visco_sl*Aw*(Pw(i+1) - Pw(i))/dx*charaP;
                    Jac(i,i+1) = Jac(i,i+1) +  dens_sl*Kw/visco_sl*Aw/dx*charaP;
                    Jac(i,i)  =  Jac(i,i) - dens_sl*Kw/visco_sl*Aw/dx*charaP + dens_dp*Kw/visco_sl*Aw*(Pw(i+1) - Pw(i))/dx*charaP;
                end
                if Pw(i-1) > Pw(i)
                    mm = i - 1;
                    visco_sl = visco_slv(mm);
                    dens_slj = dens_slv(mm);
                    dens_dpj = dens_dpv(mm);
                    RHS0(i) = RHS0(i) + dens_slj*Kw/visco_sl*Aw*(Pw(i-1) - Pw(i))/dx*charaP;
                    Jac(i,i-1) = Jac(i,i-1) +  dens_slj*Kw/visco_sl*Aw/dx*charaP;
                    Jac(i,i-1) = Jac(i,i-1) +  dens_dpj*Kw/visco_sl*Aw*(Pw(i-1) - Pw(i))/dx*charaP;
                    Jac(i,i)  =  Jac(i,i) - dens_slj*Kw/visco_sl*Aw/dx*charaP;
                    Vel(i) = Kw/visco_sl*(Pw(i-1) - Pw(i))/dx*charaP;
                else
                    mm = i;
                    visco_sl = visco_slv(mm);
                    dens_sl = dens_slv(mm);
                    dens_dp = dens_dpv(mm);
                    RHS0(i) = RHS0(i) + dens_sl*Kw/visco_sl*Aw*(Pw(i-1) - Pw(i))/dx*charaP;
                    Jac(i,i-1) = Jac(i,i-1) +  dens_sl*Kw/visco_sl*Aw/dx*charaP;
                    Jac(i,i)  =  Jac(i,i) - dens_sl*Kw/visco_sl*Aw/dx*charaP + dens_dp*Kw/visco_sl*Aw*(Pw(i-1) - Pw(i))/dx*charaP;
                    Vel(i) = Kw/visco_sl*(Pw(i-1) - Pw(i))/dx*charaP;
                end
                
            end
            % Accum Term
            dens_sl = dens_slv(i);
            dens_dp = dens_dpv(i);
            dens_sl0 = dens_slv0(i);
            RHS0(i) = RHS0(i) - dx*Aw*(dens_sl - dens_sl0)/dt;
            Jac(i,i) = Jac(i,i) - dens_dp*dx*Aw/dt;
        end
        
        dP = -Jac\RHS0;
        Pw = Pw + dP;
        epss = norm(dP)/max(Pw);
        epsr = norm(RHS0);
        eps = max(epss,epsr);
    end
    
    %MassBalance Check
    M = 0;
    for i = 1 :nG
        xfluid = Xf(i:nG:(nFluid-1)*nG+i);
        cprop =  Cp(i:nG:(nProp-1)*nG+i);
        [dens,~] = calcSLdens(Pw(i)*charaP,nFluid,fdens,fcmp,frefp,xfluid,cprop,propdens);
        xfluid = Xf0(i:nG:(nFluid-1)*nG+i);
        cprop = Cp0(i:nG:(nProp-1)*nG+i);
        [dens0,~] = calcSLdens(Pw0(i)*charaP,nFluid,fdens,fcmp,frefp,xfluid,cprop,propdens);
        M  = M + (dens - dens0)*Aw*dx;
    end
    for i = 1 : nPerf
        numPf = Perf(i);
        numFrac = well{welli}.Perfindex(numPf,:);
        if ActEle(numFrac(1)) > 0.1
            M = M+well{welli}.Sch(Schindex).PfQsl(i,1)*dt;
        end
        if ActEle(numFrac(2)) > 0.1
            M = M+well{welli}.Sch(Schindex).PfQsl(i,2)*dt;
        end
    end
    epsM = (M - mass_rate*Aw*dt)/M;
    % Update Xf and Cp by individual
    % Proppant
    MatCp =MatCp*0;
    RhsCp =RhsCp*0;
    for j = 1 : nProp
        % Flux from well to fracture
        for i = 1 : nPerf
            numPf = Perf(i);
            numFrac = well{welli}.Perfindex(numPf,:);
            numWell = well{welli}.PerfWellindex(numPf);
            if ActEle(numFrac) > 0.1
                vel = (well{welli}.Sch(Schindex).Pf_Q(i,1)/Aw+well{welli}.Sch(Schindex).Pf_Q(i,2)/Aw);
                MatCp(numWell+(j-1)*nG,numWell+(j-1)*nG) = MatCp(numWell+(j-1)*nG,numWell+(j-1)*nG) - vel/dx;
            end
        end
        for i = 1 : nG
            if i > 1
                if Vel(i) > 0
                    MatCp(i+(j-1)*nG,i-1+(j-1)*nG) =MatCp(i+(j-1)*nG,i-1+(j-1)*nG) +  Vel(i)/dx;
                else
                    MatCp(i+(j-1)*nG,i+(j-1)*nG) = MatCp(i+(j-1)*nG,i+(j-1)*nG) + Vel(i)/dx;
                end
            end
            if i < nG
                if Vel(i+1) < 0
                    MatCp(i+(j-1)*nG,i+1+(j-1)*nG) = MatCp(i+(j-1)*nG,i+1+(j-1)*nG) - Vel(i+1)/dx;
                else
                    MatCp(i+(j-1)*nG,i+(j-1)*nG) = MatCp(i+(j-1)*nG,i+(j-1)*nG) - Vel(i+1)/dx;
                end
            end
            MatCp(i+(j-1)*nG,i+(j-1)*nG) = MatCp(i+(j-1)*nG,i+(j-1)*nG)-1/dt;
            RhsCp(i+(j-1)*nG) = RhsCp(i+(j-1)*nG) - Cp0(i+(j-1)*nG)/dt;
            if j == propindex && i == 1
                RhsCp(i+(j-1)*nG) = RhsCp(i+(j-1)*nG) - rate/dx*well{welli}.Sch(Schindex).PropFraction;
            end
            if abs(MatCp(i+(j-1)*nG,i+(j-1)*nG)) < 1e-6
                MatCp(i+(j-1)*nG,i+(j-1)*nG)  = 1;
            end
        end
    end
    Cp2 = MatCp\RhsCp;
    epsCp = norm(Cp - Cp2)/sum(Cp2);
    %  Cp = Cp2;
    % M = 0;
    % for i = 1 :nG
    %     cprop =  Cp(i:nG:(nProp-1)*nG+i);
    %     cprop0  = Cp0(i:nG:(nProp-1)*nG+i);
    %     M  = M + (cprop - cprop0)*Aw*dx;
    % end
    % epsCp = (M(1) - rate*well{welli}.Sch(Schindex).PropFraction*Aw*dt)/M(1);
    %% Update Xf
    
    MatXf =MatXf*0;
    RhsXf =RhsXf*0;
    
    test  = RhsXf;
    for j = 1 : nFluid
        % Flux from well to fracture
        for i = 1 : nPerf
            numPf = Perf(i);
            numFrac = well{welli}.Perfindex(numPf);
            numWell = well{welli}.PerfWellindex(numPf);
            cprop =  Cp(numWell:nG:(nProp-1)*nG+numWell);
            densf = densfv(numWell,j);
            if ActEle(numFrac) > 0.1
                vel = (well{welli}.Sch(Schindex).Pf_Q(i,1)/Aw+well{welli}.Sch(Schindex).Pf_Q(i,2)/Aw);
                MatXf(numWell+(j-1)*nG,numWell+(j-1)*nG) =MatXf(numWell+(j-1)*nG,numWell+(j-1)*nG) - (1-sum(cprop))*densf*vel/dx;
            end
        end
        for i = 1 : nG
            if i > 1
                if Vel(i) > 0
                    mm = i-1;
                    cprop =  Cp(mm:nG:(nProp-1)*nG+mm);
                    densf = densfv(mm,j);
                    MatXf(i+(j-1)*nG,i-1+(j-1)*nG) =MatXf(i+(j-1)*nG,i-1+(j-1)*nG) + (1-sum(cprop))*densf*Vel(i)/dx;
                    test(i+(j-1)*nG) = test(i+(j-1)*nG) + (1-sum(cprop))*densf*Vel(i)/dx;
                else
                    mm = i;
                    cprop =  Cp(mm:nG:(nProp-1)*nG+mm);
                    densf = densfv(mm,j);
                    MatXf(i+(j-1)*nG,i+(j-1)*nG) = MatXf(i+(j-1)*nG,i+(j-1)*nG) + (1-sum(cprop))*densf*Vel(i)/dx;
                    test(i+(j-1)*nG) = test(i+(j-1)*nG) + (1-sum(cprop))*densf*Vel(i)/dx;
                end
            end
            if i < nG
                if Vel(i+1) < 0
                    mm = i+1;
                    cprop =  Cp(mm:nG:(nProp-1)*nG+mm);
                    densf = densfv(mm,j);
                    MatXf(i+(j-1)*nG,i+1+(j-1)*nG) = MatXf(i+(j-1)*nG,i+1+(j-1)*nG) - (1-sum(cprop))*densf*Vel(i+1)/dx;
                    test(i+(j-1)*nG) = test(i+(j-1)*nG) + (1-sum(cprop))*densf*Vel(i+1)/dx;
                else
                    mm = i;
                    densf = densfv(mm,j);
                    MatXf(i+(j-1)*nG,i+(j-1)*nG) = MatXf(i+(j-1)*nG,i+(j-1)*nG) - (1-sum(cprop))*densf*Vel(i+1)/dx;
                    test(i+(j-1)*nG) = test(i+(j-1)*nG) - (1-sum(cprop))*densf*Vel(i+1)/dx;
                end
            end
            cprop =  Cp(i:nG:(nProp-1)*nG+i);
            cprop0 =  Cp0(i:nG:(nProp-1)*nG+i);
            densf = densfv(i,j);
            densf0 = densf0v(mm,j);
            MatXf(i+(j-1)*nG,i+(j-1)*nG) = MatXf(i+(j-1)*nG,i+(j-1)*nG)-(1-sum(cprop))*densf*1/dt;
            RhsXf(i+(j-1)*nG) = RhsXf(i+(j-1)*nG) - (1-sum(cprop0))*densf0*Xf0(i+(j-1)*nG)/dt;
            test(i+(j-1)*nG) = test(i+(j-1)*nG) - ((1-sum(cprop))*(densf)*Xf(i+(j-1)*nG) - (1-sum(cprop0))*densf0*Xf0(i+(j-1)*nG))/dt;
            if j == fluidindex && i == 1
                cprop = zeros(nProp,1);
                cprop(propindex) =  well{welli}.Sch(Schindex).PropFraction;
                RhsXf(i+(j-1)*nG) = RhsXf(i+(j-1)*nG) - (1-sum(cprop))*rate/dx*dens_flinj;
                test(i+(j-1)*nG) = test(i+(j-1)*nG) + (1-sum(cprop))*rate/dx*dens_flinj;
            end
            if abs(MatXf(i+(j-1)*nG,i+(j-1)*nG)) < 1e-6
                MatXf(i+(j-1)*nG,i+(j-1)*nG)  = 1;
            end
        end
    end
    Xf2 = MatXf\RhsXf;
    epsxf = 0;
    temp = abs(max(Xf2 - Xf));
    if temp > 1e-2
        Xf = Xf2;
        epsxf = norm(Xf2 - Xf)/sum(Xf);
    end
    eps2 = max(epsxf, epsCp);
    if iIter2 > 1.1
        Pwi2 = Pw;
        Pw = (1-alphaR)*Pw + alphaR*Pwi;
        Pwi = Pwi2;
    else
        Pwi = Pw;
    end
    Cp = Cp2;
end
for i = 1 : nPerf
    for mm = 1 : nProp
        well{welli}.Sch(Schindex).PfCp(i,mm) = Cp(nG*(mm-1)+mm);
    end
    for mm = 1 : nFluid
        well{welli}.Sch(Schindex).PfXf(i,mm) = Xf(nG*(mm-1)+mm);
    end
end
% Xf = Xf2;
%MassBalance Check
M = 0;
for i = 1 :nG
    xfluid = Xf(i:nG:(nFluid-1)*nG+i);
    cprop =  Cp(i:nG:(nProp-1)*nG+i);
    [dens,~] = calcSLdens(Pw(i)*charaP,nFluid,fdens,fcmp,frefp,xfluid,cprop,propdens);
    xfluid = Xf0(i:nG:(nFluid-1)*nG+i);
    cprop = Cp0(i:nG:(nProp-1)*nG+i);
    [dens0,~] = calcSLdens(Pw0(i)*charaP,nFluid,fdens,fcmp,frefp,xfluid,cprop,propdens);
    M  = M + (dens - dens0)*Aw*dx;
end
for i = 1 : nPerf
    numPf = Perf(i);
    numFrac = well{welli}.Perfindex(numPf,:);
    if ActEle(numFrac(1)) > 0.1
        M = M+well{welli}.Sch(Schindex).PfQsl(i,1)*dt;
    end
    if ActEle(numFrac(2)) > 0.1
        M = M+well{welli}.Sch(Schindex).PfQsl(i,2)*dt;
    end
end
eps = (M - mass_rate*Aw*dt)/M;
disp('Relative Mass Error :');
fprintf('%f\n',eps);
disp('Iter Num:');
fprintf('%d %d\n',iIter1,iIter2);
end