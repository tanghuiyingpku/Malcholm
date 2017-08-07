function dt =  FluidSolidFullCouple_Scaling_well(CurT,dt,isGrow,nt)
% The boundary condition for fluid flow in no flow bounday
% In this function Shear Displacements has nothing to do with fluid flow so
% is not taken into direct calculation
% The elements with no fluid flow is not taken into calculation
%% Get global variables
global Index IndexInv nAct nAllEle_global  sigmaN_global
global TipStates EleType;
global nwell well Mat Fluid isMechActive_global;
global  MaxEle ActEle;
global DD_global PresF_global AllEle_global ConnList_global;
global CpF_global XfF_global  Fractures InitialAperture
global Hslurry_global Hbanking_global e_global Tau_global einit_global MassBL;
Acum = zeros(nAct,1);
Flow = zeros(nAct*4,3);
dt0 = dt;
CL = 0;%Fluid.spurtslop;
clear isMechActive0;
hasDecrease = 0;
Schindex = zeros(nwell,1);
for welli = 1 : nwell
    for i = 1 : well{welli}.nSch
        if CurT >= well{welli}.Sch(i).t0 && CurT <= well{welli}.Sch(i).t1
            Schindex(welli) = i;
        end
    end
end

%% Open Space for variables
restart = 1;
P_pre = zeros(nAct,1);
hslurry = zeros(nAct,Fluid.nprop);
hbanking = zeros(nAct,Fluid.nprop);
DD_pre = zeros(nAct*2,1);
Cp0 = zeros(nAct,Fluid.nprop);
Xf0 = zeros(nAct,Fluid.nfluid);
Tipstate = zeros(nAct,1);
% Assume the maximum connection number is 4
Vel_sl = zeros(nAct,5);
visco_flv = zeros(nAct,1);
visco_slv = zeros(nAct,1);
dens_slv = zeros(nAct,1);
dens_slv0 = zeros(nAct,1);
dens_dpv = zeros(nAct,1);
densfv = zeros(nAct,Fluid.nfluid);
densf0v = zeros(nAct,Fluid.nfluid);
fvisco = zeros(Fluid.nfluid,1);
fdens = zeros(Fluid.nfluid,1);
frefp = zeros(Fluid.nfluid,1);
fcmp = zeros(Fluid.nfluid,1);
propdens = zeros(Fluid.nfluid,1);
InsituS = zeros(nAct*2,1);
AllEle = zeros(nAct,11);
ConnList = zeros(nAct,8);

Q = zeros(50,1);
unitQ = 1;
unitC = 1;
epsd = 1e-8;
Sxx = Mat.Sxx*1e6;
Syy = Mat.Syy*1e6;
Sxy = Mat.Sxy*1e6;


%Scaling characteristic value
Pchara = 1e6;
Dchara = 1e-3;
Vchara = 1e-3*Mat.h*0.1;
Maxtau = 12*1e6;
einit = zeros(nAct,1);
%% Assign values
for i = 1 : nAct
    DD_pre(i) = DD_global(IndexInv(i));
    DD_pre(i+nAct) = DD_global(MaxEle+IndexInv(i));
    P_pre(i) = PresF_global(IndexInv(i))*1e6;
    Cp0(i,:) = CpF_global(IndexInv(i),:);
    Xf0(i,:) = XfF_global(IndexInv(i),:);
    Tipstate(i) = TipStates(IndexInv(i),:);
    hslurry(i,:) = Hslurry_global(IndexInv(i),:);
    hbanking(i,:) = Hbanking_global(IndexInv(i),:);
    AllEle(i,:) = AllEle_global(IndexInv(i),:);
    ConnList(i,:) = ConnList_global(IndexInv(i),:);
    %alpha = theta + 90;
    sinalp = AllEle(i,6);  %cosbet(i);
    cosalp = -AllEle(i,5); %-sinbet(i);
    InsituS(i) = (Sxx - Syy) * cosalp * sinalp - Sxy*(cosalp^2 - sinalp^2);
    InsituS(i + nAct) = Sxx * cosalp^2 + Syy * sinalp^2+2*Sxy*sinalp*cosalp;
    einit(i) = e_global(IndexInv(i));
end

hslurry0 = hslurry;
if EleType == 1
    if  isGrow > 0.1
        [extraBC,SS]  = BuildCoefMatix_Constant2(Mat,Fractures,AllEle_global(1:nAllEle_global,:),ConnList_global(1:nAllEle_global,:),1);
    else
        [extraBC,SS]  =BuildCoefMatix_Constant2(Mat,Fractures,AllEle_global(1:nAllEle_global,:),ConnList_global(1:nAllEle_global,:),0);
    end
else
    if isGrow > 0.1
        [extraBC,SS]  = BuildCoefMatix_H2(Mat,Fractures,AllEle_global(1:nAllEle_global,:),ConnList_global(1:nAllEle_global,:),1);
    else
        [extraBC,SS]  = BuildCoefMatix_H2(Mat,Fractures,AllEle_global(1:nAllEle_global,:),ConnList_global(1:nAllEle_global,:),0);
    end
end


%friction coefficient
fric = tand(Mat.fric);
[isMechActive0,ShearBC,sigmaN,e0,tau] = CheckStates(einit,nAct, AllEle, epsd,DD_pre,P_pre,Mat, fric, Maxtau);
% isMechActive0 = isMechActive0*0+1;
e = e0;
for i =  1 : nAct
    if abs(einit(i) - e(i)) > 1e-8 %%&& AllEle(i,10) > 1.1 && einit(i) < 1e-8
        einit_global(IndexInv(i)) = e(i);
    end
end
BC_insitu = InsituS;
hbanking0 = hbanking;

% Give initial opening to tip elements
DDcur = DD_pre;
nQ = 0;

% Timestep range
leftT = -1e6;
rightT = 1e6;

%value of well unkowns initials
nQ0 = 0;
for im = 1 : nwell
    wellOpen = 0;
    if Schindex(im) < 0.1 || isempty(well{im}.Sch(Schindex(im)))
        continue;
    end
    for jj = 1 :well{im}.Sch(Schindex(im)).nPf
        iele =  well{im}.Sch(Schindex(im)).Pf(jj);
        num = Index(well{im}.Perfindex(iele,:));
        if min(num) > 1e-6
            if wellOpen < 0.1
                nQ = nQ + 1;
                Q(nQ) = -Mat.Sxx*1.2*1e6/Pchara;
                wellOpen = 1;
                nQ = nQ + 1;
                Q(nQ) = well{im}.Sch(Schindex(im)).Pf_Q(jj,1)*2;
            else
                nQ = nQ + 1;
                Q(nQ) = well{im}.Sch(Schindex(im)).Pf_Q(jj,1)*2;
            end
            
        end
    end
    nQ_constv(im) = nQ-nQ0;
    nQ0 = nQ;
end
Q = Q(1:nQ);
Qinit = Q;

% For viscosity calculation
cmax = 0.7;
n = 1.3;
for i = 1 : Fluid.nfluid
    fvisco(i) = Fluid.fluid{i}.visco;
    fdens(i) = Fluid.fluid{i}.dens;
    frefp(i) = Fluid.fluid{i}.refp;
    fcmp(i) = Fluid.fluid{i}.compr;
end
%
for i = 1 : Fluid.nprop
    propdens(i) = Fluid.prop{i}.dens;
end

Khf = 1e-5;
% Convergence Criteria
tol = 1e-4;
Cp = Cp0;
Xf = Xf0;
%Previous Time Step Values
iIterT = 0;
% Initial values
Jac0 = zeros(2*nAct + nQ,2*nAct + nQ);
RHS0 = zeros(2*nAct + nQ,1);
BC0   = zeros(2*nAct,1);
ShearBC0 = ShearBC;
isMechActive = isMechActive0;

for im = 1 : nAct
    xfluid0 = Xf0(im,:);
    cprop0 =  Cp0(im,:);
    cprops = sum(cprop0);
    visco_flv(im)= dot(xfluid0,fvisco);
    visco_slv(im) = visco_flv(im)*(1-cprops/cmax)^-n;
    for mm = 1 : Fluid.nfluid
        densf0v(im,mm) = CalcFLDens(mm,P_pre(im)/Pchara);
    end
    % Original RHS
    [dens_slv0(im),~] = calcSLdens(P_pre(im),Fluid.nfluid,fdens,fcmp,frefp,xfluid0,cprop0,propdens);
end
hasIncrease = 0;
while restart > 0.1 || iIterT > 0.1
    % isMechActive = isMechActive0;
    isMechActive = isMechActive0;
    %     if iIterT > 0.1
    %         [~,~,~,isMechActive] = updateStates(isMechActive,AllEle,P*Pchara,nAct, epsd,Ds,Ds0,Dn,Mat, fric, Maxtau);
    %     end
    [~,~,sigmaN,~] = CheckStates(e0,nAct, AllEle, epsd,DD_pre,P_pre,Mat, fric, Maxtau);
    
    %n Ds n Dn n Pressure
    % Project the in-situ stresses to fractures;
    restart = 0;
    % No fluid Element
    isRegroup = 1;
    % Initial values
    P = P_pre/Pchara;
    Dn = DDcur(nAct+1:2*nAct);
    Ds = DDcur(1:nAct);
    Ds0 = DDcur(1:nAct);
    DD = DDcur;
    %while isRegroup > 0.1
    dis = 1e6;
    iIter = 0;
    % Zeros
    RHS = RHS0;
    BC = BC0;
    Jac = Jac0;
    e = e0;
    ShearBC = ShearBC0;
    % Update e
    Q = Qinit;
    Ds0i = Ds0;
    while dis > tol || isRegroup > 0.1
        Acum = Acum *0;
        Flow = Flow *0;
        % Fluid Properties
        for im = 1 : nAct
            xfluid = Xf(im,:);
            cprop =  Cp(im,:);
            for mm = 1 : Fluid.nfluid
                densfv(im,mm) = CalcFLDens(mm,P(im));
            end
            % Original RHS
            [dens_slv(im),dens_dpv(im)] = calcSLdens(P(im)*Pchara,Fluid.nfluid,fdens,fcmp,frefp,xfluid,cprop,propdens);
        end
        % For Debug USE
        dens_dpv = dens_dpv*Pchara;
        %
        % isMechActive0 = isMechActive;
        RHS = RHS*0;
        BC = BC*0;
        %Solid Part DDM
        Jac = Jac * 0 ;
        iIter = iIter + 1;
        % Pure Solid Part Solid-Solid
        disp('Building Influence Coefficient Matrix');
        
        
        disp('Coefficient Matrix Building finish');
        BC = BC + extraBC;
        % Elimination of Shear Displacements
        A110 = SS(1:nAct,1:nAct);
        A12 = SS(1:nAct,nAct+1:2*nAct);
        A21 = SS(nAct+1:2*nAct,1:nAct);
        A22 = SS(nAct+1:2*nAct,nAct+1:2*nAct);
        tempI = eye(nAct);
        EleL  = zeros(nAct,1);
        % For Pressure in Solid equation SF part Solid-Fluid
        constV = 1;
        for ii = 1 : nAct
            EleL(ii) = AllEle(ii,7);
            if isMechActive(ii) < 0.1 || Dn(ii) > -1e-16
                if isMechActive(ii) < -0.1
                    %  totally closed
                    BC(ii) = Ds0(ii)*constV;
                    BC(nAct+ii) = 0;
                    A110(ii,:) = 0 ;
                    A110(ii,ii) = 1*constV;
                    A12(ii,:) = 0 ;
                    A22(ii,:) = 0 ;
                    A22(ii,ii) = 1*constV;
                    A21(ii,:) = 0 ;
                else
                    % Sliding, Only Ds is included in calculation
                    %Shear stress equal to Mohr-colume  calculation
                    BC(nAct+ii) = 0;
                    A22(:,ii) = 0;
                    A12(:,ii) = 0;
                    A22(ii,:) = 0 ;
                    A22(ii,ii) = 1*constV;
                    A21(ii,:) = 0 ;
                    BC(ii) = BC(ii)-BC_insitu(ii) + ShearBC(ii);
                end
            else
                if isMechActive(ii) > 2.1
                    BC(ii) = Ds0(ii)*constV;
                    BC(nAct+ii) = 0;
                    A110(ii,:) = 0 ;
                    A110(ii,ii) = 1*constV;
                    A12(ii,:) = 0 ;
                    A22(ii,:) = 0 ;
                    A22(ii,ii) = 1*constV;
                    A21(ii,:) = 0 ;
                else
                    % Open Fracture
                    BC(ii) = BC(ii)-BC_insitu(ii);
                    BC(nAct+ii) =BC(nAct+ii) -BC_insitu(nAct+ii) - P(ii)*Pchara;
                end
            end
        end
        invA11 = A110\tempI;
        R1 = BC(1:nAct);
        R2 = BC(nAct+1:2*nAct);
        A11 = -A21*invA11*A12 + A22;
        RD = R2 - A21*invA11*R1;
        %Only Include the Active Grids
        %Elements with fluid
        Jac(1:nAct,1:nAct) = A11;
        RDs = RD;
        Dns = Dn;
        Pact = P;
        %Derivative of pressure
        for ii = 1 : nAct
            if isMechActive(ii) > 0.1
                %-- Scaling--
                Jac(ii,nAct+ii) = 1*Pchara;
            end
            
        end
        % Unit MPa, SI unit
        RHS(1:nAct) =A11*Dns-RDs;
        %----------Scaling------------
        Jac(1:nAct,1:nAct) = Jac(1:nAct,1:nAct) * Dchara;
        %-------------------------------
        for ii = 1 : nAct
            if isMechActive(ii) < 1e-6 || isMechActive(ii) > 2.1 %% New Type
                % sliding  or closed
                RHS(ii) = 0;
                Jac(ii,:) = 0;
                Jac(ii,ii) = 1E6;
            end
        end
        
        h = Mat.h;
        %Normal Opening :Open--> is negative
        if min(e) < 1e-16
            %  keyboard;
            ('Negative Opening');
        end
        %Volumn
        
        disp('Building Jacob Matrix of Fluid Flow Equations');
        numC = 0;
        %% Flow equations
        for i = 1: nAct
            if e(i) < 1e-16
                e(i)= 0;
            end
            if AllEle_global(IndexInv(i),10) > 1.1
                %Shearing opened element
                if isMechActive(i) < -0.1
                    [e(i),~] = calcNfWidth_S(Ds0(i), sigmaN(i),P(i)*Pchara,InsituS(i + nAct));
                end
                if isMechActive(i) < 0.1
                    [e(i),~] = calcNfWidth_S( Ds0(i),sigmaN(i),P(i)*Pchara,InsituS(i + nAct));
                end
                if isMechActive(i) > 0.1
                    [ei,~] = calcNfWidth_S(Ds0(i),0,0,InsituS(i + nAct));
                    e(i) = ei - Dn(i);
                end
            end
        end
        V = e.*EleL*Mat.h;
        for i = 1 : nAct
            h = Mat.h - sum(hbanking(i,:));
            h0 = Mat.h - sum(hbanking0(i,:));
            %Flow
            ki = max(e(i)^2/12,InitialAperture^2/12);
            nConn = ConnList(i,2);
            if AllEle_global(IndexInv(i),10) > 1.1
                % opened NF element
                if isMechActive(i) > 0.1
                    [ei,~] = calcNfWidth_S(Ds0(i),0,0,InsituS(i + nAct));
                    ki = (Dn(i))^2/12 +ei*Khf;
                else
                    % Sliding NF element
                    [ei,~] = calcNfWidth_S(Ds0(i), sigmaN(i),P(i)*Pchara,InsituS(i + nAct));
                    ki = ei*Khf;
                end
            end
            
            Li = AllEle(i,7);
            % Flow Flux Term
            
            for j = 1 : nConn
                if ConnList(i,j+2) > -0.1
                    Conj =  ConnList(i,j+2);
                    if Conj < 0.1 || Index(Conj) < 1E-9
                        continue;
                    end
                    Lj = AllEle_global(Conj,7);
                    dL = (Li + Lj)/2;
                    kj = max(e(Index(Conj))^2/12,InitialAperture^2/12);
                    % NF
                    if AllEle_global((Conj),10) > 1.1
                        % opened NF element
                        if isMechActive(Index(Conj)) > 0.1
                            [ei,~] = calcNfWidth_S(Ds0(Index(Conj)),0,0,InsituS(Index(Conj) + nAct));
                            ej = -Dn(Index(Conj))+ ei;
                            kj = Dn(Index(Conj))^2/12 +ei*Khf;
                        else
                            % Sliding NF element
                            [ei,~] = calcNfWidth_S( Ds0(Index(Conj)), sigmaN(Index(Conj)),P(Index(Conj))*Pchara,InsituS(Index(Conj) + nAct));
                            ej = ei;
                            kj = ei*Khf;
                        end
                    else
                        ej = e(Index(Conj));
                    end
                    e(Index(Conj)) = ej;
                    A = (e(i)+ej)/2*h;
                    kij = (ki*Li+kj*Lj)/(Li+Lj);
                    if isMechActive(Index(Conj)) > 0.1
                        dkj_dj = -e(Index(Conj))/6;
                        dA_dj = -h/2;
                        T_difDj = dA_dj*kij + A*(dkj_dj*Lj)/(Li+Lj);
                        T_difPj = 0;
                    else
                        T_difDj =0;
                        [~,dwdp] = calcNfWidth_S(  Ds0(Index(Conj)),sigmaN(Index(Conj)),P(Index(Conj))*Pchara,InsituS(Index(Conj) + nAct));
                        dwdp = dwdp * Pchara;
                        dkj_dPj = Khf*dwdp;
                        dA_dPj = 0;
                        T_difPj = dA_dPj*kij + A*(dkj_dPj*Lj)/(Li+Lj);
                    end
                    if isMechActive(i) > 0.1
                        dki_di = -e(i)/6;
                        dA_di = -h/2;
                        T_difDi = dA_di*kij + A*(dki_di*Li)/(Li+Lj);
                        T_difPi = 0;
                    else
                        T_difDi = 0;
                        [~,dwdp] = calcNfWidth_S(Ds0(i),sigmaN(i),P(i)*Pchara,InsituS(i + nAct));
                        dwdp = dwdp * Pchara;
                        dki_dPi = Khf*dwdp;
                        dA_dPi = 0;
                        T_difPi = dA_dPi*kij + A*(dki_dPi*Li)/(Li+Lj);
                    end
                    Tij = A*kij/dL;
                    %Upwind format
                    if P(i) >= P(Index(Conj))
                        viscosl =visco_slv(i);
                        dens_sl = dens_slv(i);
                        dens_dp = dens_dpv(i);
                        Jac(nAct+i,nAct+i) = Jac(nAct+i,nAct+i) +  Pchara*1/viscosl*(Tij * (dens_dp*(P(Index(Conj)) - P(i)) - dens_sl)*unitC + T_difPi * dens_sl*(P(Index(Conj)) - P(i)));
                        Jac(nAct+i,nAct+Index(Conj)) = Jac(nAct+i,nAct+Index(Conj)) + Pchara*1/viscosl * Tij * dens_sl*unitC + Pchara*T_difPj * dens_sl*(P(Index(Conj)) - P(i));
                        Vel_sl(i,j) = 1/viscosl * kij * (P(Index(Conj)) - P(i))/dL*Pchara;
                    else
                        viscosl =visco_slv(Index(Conj));
                        dens_sl = dens_slv(Index(Conj));
                        dens_dp = dens_dpv(Index(Conj));
                        Jac(nAct+i,nAct+i) = Jac(nAct+i,nAct+i) +  Pchara*1/viscosl * (Tij* (-dens_sl)*unitC + T_difPi * dens_sl*(P(Index(Conj)) - P(i)));
                        Jac(nAct+i,nAct+Index(Conj)) = Jac(nAct+i,nAct+Index(Conj)) +Pchara*1/viscosl*( Tij * (dens_sl+dens_dp*(P(Index(Conj)) - P(i)))*unitC + T_difPj * dens_sl*(P(Index(Conj)) - P(i)));
                        Vel_sl(i,j) = 1/viscosl * kij * (P(Index(Conj)) - P(i))/dL*Pchara;
                    end
                    if isnan(Jac(nAct+i,nAct+i))
                        %                         keyboard;
                    end
                    Jac(nAct+i,i) = Jac(nAct+i,i) +1/viscosl*T_difDi * Dchara * (-1)* dens_sl*(P(Index(Conj)) - P(i))*Pchara;
                    Jac(nAct+i,Index(Conj)) = Jac(nAct+i,Index(Conj)) +1/viscosl*T_difDj* Dchara * (-1)*dens_sl*(P(Index(Conj)) - P(i))*Pchara;
                    RHS(i + nAct) = RHS(i + nAct) + 1/viscosl*Tij * dens_sl*(P(Index(Conj)) - P(i))*Pchara;
                    numC = numC + 1;
                    Flow(numC,1) = i;
                    Flow(numC,2) = Index(Conj);
                    Flow(numC,3) = 1/viscosl*Tij * dens_sl*(P(Index(Conj)) - P(i))*Pchara;
                end
            end
            dens_sl = dens_slv(i);
            dens_sl0 = dens_slv0(i);
            dens_dp = dens_dpv(i);
            if isMechActive(i) < 0.1
                [~,dwdp] = calcNfWidth_S(Ds0(i),sigmaN(i),P(i)*Pchara,InsituS(i + nAct));
                dwdp = dwdp * Pchara;
                dDdp = 0;
            else
                dwdp = 0;
                dDdp = -1;
            end
            % Accum Inj and Leak
            % Accum
            % For Pressure Derivative
            Jac(nAct+i,nAct+i) = Jac(nAct+i,nAct+i) - (V(i)*dens_dp/dt*unitQ + (dwdp*dens_sl*h*EleL(i))/dt*unitQ)*1;
            
            if isnan(Jac(nAct+i,nAct+i))
                hasDecrease = 1;
                restart = 1;
                dt =dt /3;
                break;
            end
            % For Opening Derivative
            Jac(nAct+i,i) = Jac(nAct+i,i) - dDdp*AllEle(i,7)*h*dens_sl/dt*unitQ*Dchara;
            %Newly Modified,Accumulation Terms
            dVdt = (V(i)*dens_sl-e0(i)*dens_sl0*h0*AllEle(i,7))/dt*unitQ;
            Acum(i) = dVdt;
            if abs(dVdt) < 1e-16
                dVdt = 0;
            end
            RHS(i+nAct) = RHS(i+nAct) - dVdt;
            %Carter's Leak
            RHS(i+nAct) = RHS(i+nAct) - CL*Li*h/sqrt(CurT+dt - Tau_global(IndexInv(i)))*2*dens_sl; % Two Faces Leak off
        end
        if restart > 0.1
            break;
        end
        % Well Part
        nQ = 0;
        for ii = 1 : nwell
            wellOpen = 0;
            if Schindex(ii) < 0.1 || isempty(well{ii}.Sch(Schindex(ii)))
                continue;
            end
            for jj = 1 :well{ii}.Sch(Schindex(ii)).nPf
                iele =  well{ii}.Sch(Schindex(ii)).Pf(jj);
                num = Index(well{ii}.Perfindex(iele,:));
                if min(num) > 1e-6
                    if wellOpen < 0.1
                        %% Jacob for well
                        if ii > 1.1
                            nQ = 1 + sum(nQ_constv(1:ii-1));
                        else
                            nQ = nQ + 1;
                        end
                        dens_slurry = ConstructWell();
                        RHS(num(1)+nAct) = RHS(num(1)+nAct) + Q(nQ+1)/2*dens_slurry;
                        RHS(num(2)+nAct) = RHS(num(2)+nAct) + Q(nQ+1)/2*dens_slurry;
                        wellOpen = 1;
                    else
                        nQ = nQ + 1;
                        RHS(num(1)+nAct) = RHS(num(1)+nAct) + Q(nQ+1)/2*dens_slurry;
                        RHS(num(2)+nAct) = RHS(num(2)+nAct) + Q(nQ+1)/2*dens_slurry;
                    end
                    
                end
            end
        end
        if restart > 0.1
            break;
        end
        X0 = [Dns/Dchara;Pact;Q];
        disp('Solve Linear Equation System');
        % Jac = DumpJacob(Jac);
        dx = -Jac\RHS;
        disp('Linear Equation Solved');
        X = X0 + dx;
        Dn = X(1:nAct)*Dchara;
        Pact  = X(nAct+1:nAct+nAct)*Pchara;
        Q = X(nAct*2+1:nAct*2+nQ+1);
        P = Pact/Pchara;
        Ds = invA11*R1 - invA11*A12*Dn;
        Ds1 = Ds;
        for i = 1 : nAct
            if AllEle_global(IndexInv(i),10) > 1.1
                if sign(Ds(i)-Ds0(i))*ShearBC(i) > 1e-6 && isMechActive(i) < 0.1 && abs(Ds(i)-Ds0(i)) > 1e-8
                    isMechActive(i) = -2;
                    dis = 1e6;
                    Ds(i) = Ds0(i);
                    Dn(i) = 0;
                end
                if Dn(i) > 1e-20 && isMechActive(i) > 0.1
                    isMechActive(i) = 0;
                    Dn(i) = 0;
                end
                
            else
                if Dn(i) > 1e-20
                    Dn(i) = 0;
                    %isMechActive(i) = 3;
                end
            end
        end
        Ds0i = Ds;
        e1 = e;
        [e,sigmaN,ShearBC,~] = updateStates(isMechActive,AllEle,P*Pchara,nAct, epsd,Ds1,Ds0,Dn,Mat, fric, Maxtau);
        disp('   ');
        disD = norm(RHS(1:nAct));
        disP = norm(RHS(nAct+1:nAct+nAct))/1e1;
        disQ = norm(RHS(nAct*2+1:nAct*2+nQ+1));
        xnormD = norm(dx(1:nAct))/(max(abs(X(1:nAct))));
        xnormP = norm(dx(nAct:nAct+nAct))/(max(abs(X(nAct:nAct+nAct))));
        xnormQ = norm(dx(nAct*2+1:nAct*2+nQ+1))/(max(abs(X(nAct*2+1:nAct*2+nQ+1))));
        % disStates = norm(isMechActive - isMechActive0);
        %dis = max([disD,disP,xnormD,xnormP,xnormQ,disQ]);
        dis0 = dis;
        dis = max([disD/100,disP,xnormD  xnormP xnormQ disQ]);
        %         if abs(dis - dis0)/dis < 0.1 && dis > 1e-3 && haschange < 0.1
        %             [~,~,~,isMechActive] = updateStates(isMechActive,AllEle,P*Pchara,nAct, epsd,Ds,Ds0,Dn,Mat, fric, Maxtau);
        %             haschange = 1;
        %         end
        if abs(dis - dis0)/dis < 1e-5
            dis = 0;
        end
        DD = [Ds;Dn];
        fprintf('**Number of Iteration is : %d  Relative Convergence Norm Value times %f \n', iIter,dis);
        isRegroup = 0;
        disp('Check Elements Types...');
        disp('   ');
        fprintf('If Regroup ?: %d\n',isRegroup);
        if iIter >= 30 && min(e) > 1e-20
            dt = dt/3;
            hasDecrease = 1;
            restart = 1;
            break;
        end
        if max(P) > 90
            dt = dt/3;
            hasDecrease = 1;
            restart = 1;
            break;
        end
        if iIter >= 30 && min(e) < 1e-20
            dt = dt * 1.5;
            restart = 1;
            break;
        end
        if abs(dis) < 1e-20 && hasDecrease < 0.1
%             dt = dt * 1.5;
%             if dt > 2
%                 restart = 0;
%             else
%                 restart = 1;
%             end
            %hasIncrease = 1;
            break;
        end
    end
    %% Check Time step
    temp = dt;
    if  restart < 0.1
        [KI1,KI2] = StressIntensF_local(DD);
        [~,CritK] = FindCritDir_Anisotropy(KI1,KI2);
        temp = dt;
        vmax= max(max(abs(Vel_sl)));
        pmax = max(dx(nAct+1:nAct+nAct));
        [leftT,rightT,isAdjust,dt]  = AdjustPropSpeed(leftT,rightT,dt,CritK,KI1,KI2,vmax,pmax);
        if isAdjust > 0.1
            iIterT = iIterT + 1;
        else
            iIterT = 0;
            dt = temp;
        end
        if hasDecrease > 0.1 || dt > 2
            iIterT = 0;
            dt = temp;
        end
    end
    if iIterT > 2
        iIterT = 0;
        dt = temp;
    end
    
end
isMechActive_global(1:nAct) = isMechActive;
%% check mass balance
e = e1;
detQ = MassBalance(CurT,e0,e,Schindex,dens_slv0,dens_slv,h,dt,CL);
if abs(detQ) > 1e-4
    f = 1;
    MassBL(nt) = detQ;
end
MassBL(nt) = detQ;
fprintf('**************\n');
fprintf('Relative Mass Balance = : %f\n',detQ);
% end
disp('Current Time step Coupling Equation Converged');
% Give Value to global  Variables
[Vprop_z,Vprop_x] = CalcProppantV(nAct,Cp,Xf,densfv,propdens,visco_slv,e,Vel_sl);
Cp_new = UpdateCp_x(Cp0,nAct,Fluid.nprop,Mat.h,hbanking,hbanking0,ConnList,Vprop_x,Schindex,e,e0,dt,AllEle(:,7));
Cp = Cp_new;
nprop = Fluid.nprop;
for mm = 1 : nAct
    bank_height = 0;
    bank_height0 = 0;
    for i = 1 : nprop
        num = (i-1)*nprop+mm;
        if Cp_new(mm,i) > 1e-9
            hslurry(mm,i) = hslurry0(mm,i) ;%- Vprop_z(num)*dt;
        end
        hbanking(mm,i) = hbanking0(mm,i) + Vprop_z(num)*dt*Cp_new(mm,i);
        bank_height = bank_height + hbanking(mm,i);
        bank_height0 = bank_height0 + hbanking0(mm,i);
    end
    for i = 1 : nprop
        if Cp(mm,i) > 1e-16
            Cp_new(mm,i) = (Cp(mm,i)*(h - bank_height0)*e(mm) - (hbanking(mm,i)-hbanking0(mm,i))*e(mm))/e(mm)/(h-bank_height);
        end
    end
end
Cp = Cp_new;
%% give value back to global variabls
for i = 1 : nAct
    PresF_global(IndexInv(i)) = P(i)*Pchara/1e6;
    if Dn(i) > -1e-14
        Dn(i) = 0;
    end
    if abs(Ds(i)) < 1e-20
        Ds(i) = Ds0(i)*0.5*1E-12;
    end
    if isMechActive(i) < 0.1
        DD_global(IndexInv(i)) = Ds(i);
        DD_global(IndexInv(i)+MaxEle)= 0;
    else
        DD_global(IndexInv(i))=  Ds(i);
        DD_global(IndexInv(i) + MaxEle) = Dn(i);
    end
    
    CpF_global(IndexInv(i),:) = Cp(i,:);
    XfF_global(IndexInv(i),:) = Xf(i,:);
    Hslurry_global(IndexInv(i),:) = hslurry(i,:);
    Hbanking_global(IndexInv(i),:) = hbanking(i,:);
    sigmaN_global(IndexInv(i),:) = sigmaN(i);
    e_global(IndexInv(i)) =e(i);
end
%% Give Value to well
nQ = 0;
for ii = 1 : nwell
    wellOpen = 0;
    for jj = 1 :well{ii}.Sch(Schindex(ii)).nPf
        iele =  well{ii}.Sch(Schindex(ii)).Pf(jj);
        num = Index(well{ii}.Perfindex(iele,:));
        if min(num) > 1e-6
            if wellOpen < 0.1
                dens_slurry = WellSlurryDens();
                nQ = nQ + 1;
                nQ = nQ + 1;
                well{ii}.Sch(Schindex(ii)).Pf_Q(jj,1) = Q(nQ)/2;
                well{ii}.Sch(Schindex(ii)).Pf_Q(jj,2) = Q(nQ)/2;
                well{ii}.Sch(Schindex(ii)).PfQsl(jj,1) = Q(nQ)/2*dens_slurry;
                well{ii}.Sch(Schindex(ii)).PfQsl(jj,2) = Q(nQ)/2*dens_slurry;
                well{ii}.Sch(Schindex(ii)).Pres(jj) = P(num(1));
                wellOpen = 1;
            else
                nQ = nQ + 1;
                well{ii}.Sch(Schindex(ii)).Pf_Q(jj,1) = Q(nQ)/2;
                well{ii}.Sch(Schindex(ii)).Pf_Q(jj,2) = Q(nQ)/2;
                well{ii}.Sch(Schindex(ii)).Pres(jj) = P(num(1));
                well{ii}.Sch(Schindex(ii)).PfQsl(jj,1) = Q(nQ)/2*dens_slurry;
                well{ii}.Sch(Schindex(ii)).PfQsl(jj,2) = Q(nQ)/2*dens_slurry;
            end
            
        end
    end
end
%% Well functions
    function dens_slwell = ConstructWell()
        wellD = 0.1; % Well Diameter
        nFluid = Fluid.nfluid;
        nProp = Fluid.nprop;
        proptype =  well{ii}.Sch(Schindex(ii)).Prop;
        fluidtype = well{ii}.Sch(Schindex(ii)).Fluid;
        nPerf  = well{ii}.Sch(Schindex(ii)).nPf;
        Perf = well{ii}.Sch(Schindex(ii)).Pf;
        for ifluid = 1 : nFluid
            name = Fluid.fluid{ifluid}.name;
            if strcmp(name,fluidtype) == 1
                fluidindex = ifluid;
            end
        end
        %
        for iprop = 1 : nProp
            name = Fluid.prop{iprop}.name;
            if strcmp(name,proptype) == 1
                propindex = iprop;
            end
        end
        
        xfluid_well = zeros(nFluid,1);
        xfluid_well(fluidindex) = 1;
        cprop_well = zeros(nProp,1);
        cprop_well(propindex) =  well{ii}.Sch(Schindex(ii)).PropFraction;
        isInj = 0;
        if strcmp(well{ii}.Sch(Schindex(ii)).wellType,'INJ') == 1
            isInj = 1;
        end
        if strcmp(well{ii}.Sch(Schindex(ii)).wellType,'PROD') == 1
            isInj = -1;
        end
        Q0 = isInj * well{ii}.Sch(Schindex(ii)).ContrValue*0.159/60;
        [dens_slwell,~] = calcSLdens(Mat.Pp*1e6,nFluid,fdens,fcmp,frefp,xfluid_well,cprop_well,propdens);
        Dens = dens_slwell;
        nVari = 0;
        Vari = zeros(nPerf,1);
        cpropsw = sum(cprop_well);
        visco_flvw= dot(xfluid_well,fvisco);
        visco_slvw = visco_flvw*(1-cpropsw/cmax)^-n;
        nVari = nVari + 1;
        for iperf = 1 : nPerf
            numPf = Perf(iperf);
            numFrac = well{ii}.Perfindex(numPf,:);
            if ActEle(numFrac(1)) > 0.1
                nVari = nVari + 1;
                Vari(iperf) = iperf;
            end
        end
        
        Qchara = 1;%Q0;
        Pcharaw = 1;%abs(Mat.Sxx)*1e6;
        Kd = 0.8;
        Dperf = 10*1e-3*100/2.54;
        Dens = Dens*0.0083454045;
        n2 = 6;
        % unit dens: lbm/gal dp:in RATE bpm = 0.159/60 m^3/s
        Qunit = (1*60/0.159)^2;
        K1 = 0.2369*Dens/n2^2/Dperf^4/Kd^2*Qchara^2/Pcharaw*Qunit;
        % From psi to Pa
        K1 =  K1 / 145 * 1E6;
        K1v = K1 + zeros(nVari-1,1);
        K2 = 0;%128*visco_slvw/pi/wellD^4*Qchara/Pcharaw*1;
        
        L = zeros(nVari-1,1);
        RHS(2*nAct+nQ) = sum(Q(nQ+1:nQ+nVari-1)) - Q0;
        Jac(2*nAct+nQ,2*nAct+nQ+1:2*nAct+nQ+nQ_constv(ii)-1) = 1;
        
        for ip = 1 : nVari-1
            K1i = K1v(ip);
            num_perf = Index(well{ii}.Perfindex(Perf(ip),:));
            RHS(2*nAct+nQ+ip) = Q(nQ)*Pchara -  P(num_perf(1))/Pcharaw*Pchara - K1i*Q(nQ+ip)^2;
            %Derivative of Pressure in fracture
            Jac(2*nAct+nQ +ip,nAct + num_perf(1)) =  -  1/Pcharaw*Pchara ;
            
            %Derivative of Rate in fracture flow equation
            Jac(nAct + num_perf(1),2*nAct+nQ +ip)   =  0.5 * dens_slwell;
            Jac(nAct + num_perf(2),2*nAct+nQ +ip)   =  0.5 * dens_slwell;
            
            % Derivative of Rate in well equation
            Jac(2*nAct+nQ +ip,2*nAct+nQ+ip) = Jac(2*nAct+nQ +ip,2*nAct+nQ +ip) - 2*K1i*Q(nQ+ip);
            % Derivative of well head pressure in well equation
            Jac(2*nAct+nQ +ip,2*nAct+nQ) = Pchara;
            
            if ip == 1
                L(ip) = CalculateDis(well{ii}.heel', well{ii}.Perf(ip,:));
                dPfric = K2*L(ip)*Q0/Qchara;
            else
                L(ip)  = CalculateDis(well{ii}.Perf(ip,:), well{ii}.Perf(ip-1,:));
                dPfric = K2*L(1)*Q0/Qchara;
                for jP = 1: ip-1
                    dPfric = dPfric + L(jP+1)*(Q0/Qchara - sum(Q(nQ+1:nQ+jP)))* K2;
                    Jac(2*nAct+nQ +ip, 2*nAct+nQ+jP) = Jac(2*nAct+nQ +ip, 2*nAct+nQ+jP)+ sum(L(jP+1:ip)) * K2;
                end
            end
            RHS(2*nAct+nQ+ip) =  RHS(2*nAct+nQ+ip) - dPfric;
        end
        for iq = 1 : nVari-1
            if Q(nQ+iq) < 1e-10
                % no flow
                RHS(2*nAct+nQ+iq)=0;
                Jac(2*nAct+nQ+iq,:) = 0;
                Jac(2*nAct+nQ+iq,2*nAct+nQ+iq) = 1;
                Q(nQ+iq)  =0;
            end
        end
    end
    function dens_slurry = WellSlurryDens()
        nFluid = Fluid.nfluid;
        nProp = Fluid.nprop;
        proptype =  well{ii}.Sch(Schindex(ii)).Prop;
        fluidtype = well{ii}.Sch(Schindex(ii)).Fluid;
        for ifluid = 1 : nFluid
            name = Fluid.fluid{ifluid}.name;
            if strcmp(name,fluidtype) == 1
                fluidindex = ifluid;
            end
        end
        %
        for iprop = 1 : nProp
            name = Fluid.prop{iprop}.name;
            if strcmp(name,proptype) == 1
                propindex = iprop;
            end
        end
        
        xfluid_well = zeros(nFluid,1);
        xfluid_well(fluidindex) = 1;
        cprop_well = zeros(nProp,1);
        cprop_well(propindex) =  well{ii}.Sch(Schindex(ii)).PropFraction;
        [dens_slurry,~] = calcSLdens(Mat.Pp*1e6,nFluid,fdens,fcmp,frefp,xfluid_well,cprop_well,propdens);
        
    end
%% clear temp variabls
clear e0 e P_pre hslurry hbanking DD_pre DD DDcur;
clear Cp0 Xf0 Tipstate Vel_sl visco_flv  visco_slv  dens_slv;
clear  dens_slv0  dens_dpv   densfv    densf0v
clear  fvisco fdens frefp fcmp propdens;
clear  InsituS   AllEle  ConnList sigmaN ShearBC Q BC_insitu;
clear RHS
clear Jac
end


