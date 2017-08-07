function dt =  FluidSolidFullCouple_Scaling_well1(CurT,dt,isGrow)
% The boundary condition for fluid flow in no flow bounday
% In this function Shear Displacements has nothing to do with fluid flow so
% is not taken into direct calculation
% The elements with no fluid flow is not taken into calculation
%% Get global variables
global Index IndexInv nAct nAllEle_global
global isMechActive;
global TipStates EleType;
global nwell well Mat Fluid;
global  MaxEle ActEle;
global DD_global PresF_global AllEle_global ConnList_global;
global CpF_global XfF_global  Fractures InitialAperture
global Hslurry_global Hbanking_global e_global Tau_global;
CL = Fluid.spurtslop;
clear isMechActive0;
Schindex = zeros(nwell,1);
for welli = 1 : nwell
    for i = 1 : well{welli}.nSch
        if CurT >= well{welli}.Sch(i).t0 && CurT <= well{welli}.Sch(i).t1
            Schindex(welli) = i;
        end
    end
end
restart = 1;
e0 = zeros(nAct,1);
P_pre = zeros(nAct,1);
hslurry = zeros(nAct,Fluid.nprop);
hbanking = zeros(nAct,Fluid.nprop);
DD_pre = zeros(nAct*2,1);
Cp0 = zeros(nAct,Fluid.nprop);
Xf0 = zeros(nAct,Fluid.nfluid);
Tipstate = zeros(nAct,1);
% Assume the maximum connection number is 4
Vel_sl = zeros(nAct,4);
visco_flv = zeros(nAct,1);
visco_slv = zeros(nAct,1);
dens_slv = zeros(nAct,1);
dens_slv0 = zeros(nAct,1);
dens_dpv = zeros(nAct,1);
densfv = zeros(nAct,Fluid.nfluid);
densf0v = zeros(nAct,Fluid.nfluid);

for i = 1 : nAct
    e0(i) = -DD_global(MaxEle+IndexInv(i));
    DD_pre(i) = DD_global(IndexInv(i));
    DD_pre(i+nAct) = DD_global(MaxEle+IndexInv(i));
    P_pre(i) = PresF_global(IndexInv(i))*1e6;
    Cp0(i,:) = CpF_global(IndexInv(i),:);
    Xf0(i,:) = XfF_global(IndexInv(i),:);
    Tipstate(i) = TipStates(IndexInv(i),:);
    hslurry(i,:) = Hslurry_global(IndexInv(i),:);
    hbanking(i,:) = Hbanking_global(IndexInv(i),:);
end
hbanking0 = hbanking;
% For viscosity calculation
cmax = 0.7;
n = 1.3;
fvisco = zeros(Fluid.nfluid,1);
fdens = zeros(Fluid.nfluid,1);
frefp = zeros(Fluid.nfluid,1);
fcmp = zeros(Fluid.nfluid,1);
propdens = zeros(Fluid.nfluid,1);
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
unitQ = 1;
unitC = 1;
epsd = 1e-3;
Sxx = Mat.Sxx*1e6;
Syy = Mat.Syy*1e6;
Sxy = Mat.Sxy*1e6;
%Scaling characteristic value
Pchara = 1e6;
Dchara = 1e-3;
Vchara = 1e-3*Mat.h*0.1;
%
InsituS = zeros(nAct*2,1);
AllEle = zeros(nAct,11);
ConnList = zeros(nAct,8);
for i = 1 : nAct
    AllEle(i,:) = AllEle_global(IndexInv(i),:);
    ConnList(i,:) = ConnList_global(IndexInv(i),:);
    %alpha = theta + 90;
    sinalp = AllEle(i,6);  %cosbet(i);
    cosalp = -AllEle(i,5); %-sinbet(i);
    InsituS(i) = (Sxx - Syy) * cosalp * sinalp - Sxy*(cosalp^2 - sinalp^2);
    InsituS(i + nAct) = Sxx * cosalp^2 + Syy * sinalp^2+2*Sxy*sinalp*cosalp;
end
BC_insitu = InsituS;
% Convergence Criteria
tol = 1e-4;
sigmaN = zeros(nAct,1);
ShearBC = zeros(nAct,1);
clear isMechActive0;
Cp = Cp0;
Xf = Xf0;
nQ = 0;
Q = zeros(50,1);
% Timestep range
leftT = -1e6;
rightT = 1e6;
%% Solve equations
%Size of well unkowns
for ii = 1 : nwell
    wellOpen = 0;
    for jj = 1 :well{ii}.Sch(Schindex(ii)).nPf
        iele =  well{ii}.Sch(Schindex(ii)).Pf(jj);
        num = Index(well{ii}.Perfindex(iele,:));
        if min(num) > 1e-6
            if wellOpen < 0.1
                nQ = nQ + 1;
                Q(nQ) = -Mat.Sxx*1.2*1e6/Pchara;
                wellOpen = 1;
                nQ = nQ + 1;
                Q(nQ) = well{ii}.Sch(Schindex(ii)).Pf_Q(jj,1)*2;
            else
                nQ = nQ + 1;
                Q(nQ) = well{ii}.Sch(Schindex(ii)).Pf_Q(jj,1)*2;
            end
            
        end
    end
end
Q = Q(1:nQ);
nQ_const = nQ;
% Give initial opening to tip elements
DDcur = DD_pre;
for ii = 1 : nAct
    %if TipStates(IndexInv(ii)) > 0.1
    if -(DDcur(ii+nAct) ) < InitialAperture
        DDcur(ii+nAct) = -InitialAperture;
    end
    if e0(ii) < 1e-16
        e0(ii) = InitialAperture;
    end
    %end
end
%
for ii = 1 : nAct
    isMechActive(ii) = 1;
end
isMechActive0 = isMechActive;
%Previous Time Step Values
iIterT = 0;

while restart > 0.1 || iIterT > 0.1
    %n Ds n Dn n Pressure
    isMechActive = isMechActive0;
    Jac0 = zeros(2*nAct + nQ,2*nAct + nQ);
    RHS = zeros(2*nAct + nQ,1);
    % RHS for solid equation set
    BC0   = zeros(2*nAct,1);
    % Project the in-situ stresses to fractures;
    restart = 0;
    DD = DDcur;
    P = P_pre/Pchara;
    % No fluid Element
    isRegroup = 1;
    % Regroup is for types changes
    Dn = DD(nAct+1:2*nAct);
    %Ds = DD(1:nAct);
    %while isRegroup > 0.1
    dis = 1e6;
    iIter = 0;
    BC = BC0;
    Jac = Jac0;
    for ii = 1 : nAct
        if  AllEle(ii,10) > 1.1
            % Closed Fracture
            xx= AllEle(ii,8)+epsd*AllEle(ii,5);
            yy = AllEle(ii,9)-epsd*AllEle(ii,6);
            Sxxi = 0;
            Syyi = 0;
            Sxyi = 0;
            [Sxxi,Syyi,Sxyi,~,~] = CalcPointStress_C(Mat,DD,xx,yy,AllEle_global(1:nAllEle_global,:),ConnList_global(1:nAllEle_global,:));
            Sxxi = Sxxi + Mat.Sxx;
            Syyi = Syyi + Mat.Syy;
            Sxyi = Sxyi + Mat.Sxy;
            cosb = AllEle(ii,6);
            sinb = AllEle(ii,5);
            sigmaN(ii) = abs(Sxxi*sinb^2-2*Sxyi*sinb*cosb + Syyi*cosb^2)*1e6;
            fric = tand(1);%1e-2;
            if  abs(isMechActive(ii)) < 1e-5
                ShearBC(ii) = fric*(sigmaN(ii) - P(ii)*Pchara);
            end
        end
    end
    while dis > tol || isRegroup > 0.1
        % Properties
        for ii = 1 : nAct
            xfluid = Xf(ii,:);
            xfluid0 = Xf0(ii,:);
            cprop =  Cp(ii,:);
            cprop0 =  Cp0(ii,:);
            cprops = sum(cprop);
            visco_flv(ii)= dot(xfluid,fvisco);
            visco_slv(ii) = visco_flv(ii)*(1-cprops/cmax)^-n;
            for mm = 1 : Fluid.nfluid
                densf0v(ii,mm) = CalcFLDens(mm,P(ii));
                densfv(ii,mm) = CalcFLDens(mm,P(ii));
            end
            % Original RHS
            [dens_slv(ii),dens_dpv(ii)] = calcSLdens(P(ii)*Pchara,Fluid.nfluid,fdens,fcmp,frefp,xfluid,cprop,propdens);
            [dens_slv0(ii),~] = calcSLdens(P_pre(ii),Fluid.nfluid,fdens,fcmp,frefp,xfluid0,cprop0,propdens);
        end
        % For Debug USE
        dens_dpv = dens_dpv*Pchara;
        %
        isMechActive0 = isMechActive;
        RHS = RHS*0;
        BC = BC*0;
        %Solid Part DDM
        Jac = Jac * 0 ;
        iIter = iIter + 1;
        % Pure Solid Part Solid-Solid
        disp('Building Influence Coefficient Matrix');
        if EleType == 1
            if iIter < 1.1 && isGrow > 0.1
                [extraBC,SS]  = BuildCoefMatix_Constant2(Mat,Fractures,AllEle_global(1:nAllEle_global,:),ConnList_global(1:nAllEle_global,:),1);
            else
                if iIter < 1.1
                    [extraBC,SS]  =BuildCoefMatix_Constant2(Mat,Fractures,AllEle_global(1:nAllEle_global,:),ConnList_global(1:nAllEle_global,:),0);
                end
            end
        else
            if iIter < 1.1 && isGrow > 0.1
                [extraBC,SS]  = BuildCoefMatix_H2(Mat,Fractures,AllEle_global(1:nAllEle_global,:),ConnList_global(1:nAllEle_global,:),1);
            else
                if iIter < 1.1
                    [extraBC,SS]  = BuildCoefMatix_H2(Mat,Fractures,AllEle_global(1:nAllEle_global,:),ConnList_global(1:nAllEle_global,:),0);
                end
            end
        end
        
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
        e = -Dn;
        for ii = 1 : nAct
            EleL(ii) = AllEle(ii,7);
            if isMechActive(ii) < 0.1
                if abs(isMechActive(ii)+3) < 1e-6
                    %  totally closed
                    BC(ii) = 0;
                    BC(nAct+ii) = 0;
                    A110(ii,:) = 0 ;
                    A110(ii,ii) = 1;
                    A12(ii,:) = 0 ;
                    A22(ii,:) = 0 ;
                    A22(ii,ii) = 1;
                    A21(ii,:) = 0 ;
                end
                if abs(isMechActive(ii)) < 1e-6
                    % Sliding, Only Ds is included in calculation
                    % Shear stress equal to Mohr-colume  calculation
                    BC(nAct+ii) = 0;
                    A22(ii,:) = 0 ;
                    A21(ii,:) = 0 ;
                    A12(:,ii) = 0;
                    A22(:,ii) = 0;
                    A22(ii,ii) = 1;
                    
                    ShearBC(ii) = fric*(sigmaN(ii) - P(ii)*Pchara);
                    BC(ii) = BC(ii)-BC_insitu(ii) - ShearBC(ii);
                    
                    %       BC(nAct+ii) =BC(nAct+ii) -BC_insitu(nAct+ii) - P(ii)*Pchara;
                end
            else
                BC(ii) = BC(ii)-BC_insitu(ii);
                BC(nAct+ii) =BC(nAct+ii) -BC_insitu(nAct+ii) - P(ii)*Pchara;
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
        h = Mat.h;
        %Normal Opening :Open--> is negative
        if min(e) < 1e-16
            %  keyboard;
            ('Negative Opening');
        end
        %Volumn
        V = e.*EleL*Mat.h;
        disp('Building Jacob Matrix of Fluid Flow Equations');
        for i = 1: nAct
            if e(i) < 1e-16
                e(i)= 0;
            end
            if AllEle_global(IndexInv(i),10) > 1.1
                %Shearing opened element
                if abs(isMechActive(i)+2) < 1e-6 ||  abs(isMechActive(i)+3) < 1e-6
                    [e(i),~] = calcNfWidth( sigmaN(i),P(i)*Pchara);
                end
                e(i) = max(e(i),InitialAperture);
            end
            h = Mat.h - sum(hbanking(i,:));
            h0 = Mat.h - sum(hbanking0(i,:));
            V(i) = e(i).*EleL(i)*h;
            %Flow
            nConn = ConnList(i,2);
            ki = max(e(i)^2/12,InitialAperture^2/12);
            Li = AllEle(i,7);
            % Flow Flux Term
            for j = 1 : nConn
                if ConnList(i,j+2) > -0.1
                    Conj =  ConnList(i,j+2);
                    if Conj < 0.1 || Index(Conj) < 1E-9
                        continue;
                    end
                    Lj = AllEle_global(Conj,7);
                    A = (e(i)+e(Index(Conj)))/2*h;
                    dL = (Li + Lj)/2;
                    kj = max(e(Index(Conj))^2/12,InitialAperture^2/12);
                    kij = (ki*Li+kj*Lj)/(Li+Lj);
                    if abs(isMechActive(i)+2) < 1e-6 ||  abs(isMechActive(i)+3) < 1e-6
                        T_difDj =0;
                        [~,dwdp] = calcNfWidth( sigmaN(Index(Conj)),P(Index(Conj))*Pchara);
                        dwdp = dwdp * Pchara;
                        dkj_dPj = e(Index(Conj))/6*dwdp;
                        dA_dPj = dwdp*h/2;
                        T_difPj = dA_dPj*kij + A*(dkj_dPj*Lj)/(Li+Lj);
                    else
                        dkj_dj = -e(Index(Conj))/6;
                        dA_dj = -h/2;
                        T_difDj = dA_dj*kij + A*(dkj_dj*Lj)/(Li+Lj);
                        T_difPj = 0;
                    end
                    if abs(isMechActive(i)+2) < 1e-6 ||  abs(isMechActive(i)+3) < 1e-6
                        T_difDi = 0;
                        [~,dwdp] = calcNfWidth(sigmaN(i),P(i)*Pchara);
                        dwdp = dwdp * Pchara;
                        dki_dPi = e(i)/6*dwdp;
                        dA_dPi = dwdp*h/2;
                        T_difPi = dA_dPi*kij + A*(dki_dPi*Li)/(Li+Lj);
                    else
                        dki_di = -e(i)/6;
                        dA_di = -h/2;
                        T_difDi = dA_di*kij + A*(dki_di*Li)/(Li+Lj);
                        T_difPi = 0;
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
                        keyboard;
                    end
                    Jac(nAct+i,i) = Jac(nAct+i,i) +1/viscosl*T_difDi * Dchara * (-1)* dens_sl*(P(Index(Conj)) - P(i))*Pchara;
                    Jac(nAct+i,Index(Conj)) = Jac(nAct+i,Index(Conj)) +1/viscosl*T_difDj* Dchara * (-1)*dens_sl*(P(Index(Conj)) - P(i))*Pchara;
                    RHS(i + nAct) = RHS(i + nAct) + 1/viscosl*Tij * dens_sl*(P(Index(Conj)) - P(i))*Pchara;
                end
            end
            dens_sl = dens_slv(i);
            dens_sl0 = dens_slv0(i);
            dens_dp = dens_dpv(i);
            
            if abs(isMechActive(i)+2) < 1e-6 ||  abs(isMechActive(i)+3) < 1e-6
                [~,dwdp] = calcNfWidth(sigmaN(i),P(i)*Pchara);
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
                dt = dt/3;
                restart = 1;
                break;
            end
            % For Opening Derivative
            Jac(nAct+i,i) = Jac(nAct+i,i) - dDdp*AllEle(i,7)*h*dens_sl/dt*unitQ*Dchara;
            %Newly Modified,Accumulation Terms
            dVdt = (V(i)*dens_sl-e0(i)*dens_sl0*h0*AllEle(i,7))/dt*unitQ;
            if abs(dVdt) < 1e-16
                dVdt = 0;
            end
            RHS(i+nAct) = RHS(i+nAct) - dVdt;
            %Carter's Leak
            RHS(i+nAct) = RHS(i+nAct) - CL*Li*h/sqrt(CurT+dt - Tau_global(IndexInv(i)))*2*dens_sl; % Two Faces Leak off
        end
        % Well Part
        nQ = 0;
        for ii = 1 : nwell
            wellOpen = 0;
            for jj = 1 :well{ii}.Sch(Schindex(ii)).nPf
                iele =  well{ii}.Sch(Schindex(ii)).Pf(jj);
                num = Index(well{ii}.Perfindex(iele,:));
                if min(num) > 1e-6
                    if wellOpen < 0.1
                        %% Jacob for well
                        nQ = nQ + 1;
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
        for i = 1 : nAct
            if Dn(i) > -1e-7
                Dn(i) = 0;
                % For natural fracture element
%                 if AllEle_global(IndexInv(i),10) > 1.1
%                     isMechActive(i) = 0;
%                 end
            end
        end
        DD = [Ds;Dn];
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
        dis = max([xnormD  xnormP xnormQ disQ]);
        if abs(dis - dis0) < 1e-5
            dis = 0;
        end
        fprintf('**Number of Iteration is : %d  Relative Convergence Norm Value times %f \n', iIter,dis);
        isRegroup = 0;
        disp('Check Elements Types...');
        disp('   ');
        fprintf('If Regroup ?: %d\n',isRegroup);
        %        if iIter > 5
        %       end
        if iIter >= 30
            dt = dt/3;
            restart =  1;
            break;
        end
    end
    %     for ii = 1 : nAct
    %         if  AllEle(ii,10) > 1.1 && abs(isMechActive(ii)) < 1e-5
    %             % Closed Fracture
    %             xx= AllEle(ii,8)+epsd*AllEle(ii,5);
    %             yy = AllEle(ii,9)-epsd*AllEle(ii,6);
    %             Sxxi = 0;
    %             Syyi = 0;
    %             Sxyi = 0;
    %             %[Sxxi,Syyi,Sxyi,~,~] = CalcPointStress_C(Mat,DD,xx,yy,AllEle_global(1:nAllEle_global,:),ConnList_global(1:nAllEle_global,:));
    %             Sxxi = Sxxi + Mat.Sxx;
    %             Syyi = Syyi + Mat.Syy;
    %             Sxyi = Sxyi + Mat.Sxy;
    %             cosb = AllEle(ii,6);
    %             sinb = AllEle(ii,5);
    %             sigmaN(ii) = abs(Sxxi*sinb^2-2*Sxyi*sinb*cosb + Syyi*cosb^2)*1e6;
    %             fric = tand(12);
    %             ShearBC(ii) = fric*(sigmaN(ii) - P(ii)*Pchara);
    %         end
    %     end
    %% Check Time step
    if restart < 0.1 %&& dis > 1e-10
        [KI1,KI2] = StressIntensF_local(DD);
        [~,CritK] = FindCritDir_Anisotropy(KI1,KI2);
        [leftT,rightT,isAdjust,dt]  = AdjustPropSpeed(leftT,rightT,dt,CritK,KI1,KI2);
        if isAdjust > 0.1
            iIterT = iIterT + 1;
        else
            iIterT = 0;
        end
    end
    if iIterT > 4
        iIterT = 0;
    end
end
%% check mass balance
detQ = MassBalance(CurT,e0,e,Schindex,dens_slv0,dens_slv,h,dt,CL);
fprintf('**************\n');
fprintf('Relative Mass Balance = : %f\n',detQ);
% end
disp('Current Time step Coupling Equation Converged');
% Give Value to global  Variables

%% give value back to global variabls
for i = 1 : nAct
    PresF_global(IndexInv(i)) = P(i)*Pchara/1e6;
    if Dn(i) > -1e-7
        Dn(i) = 0;
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
        Dperf = 20*1e-3*100/2.54;
        Dens = Dens*0.0083454045;
        n2 = 15;
        % unit dens: lbm/gal dp:in RATE bpm = 0.159/60 m^3/s
        Qunit = (1*60/0.159)^2;
        K1 = 0.2369*Dens/n2^2/Dperf^4/Kd^2*Qchara^2/Pcharaw*Qunit;
        % From psi to Pa
        K1 = K1 / 145 * 1E6;
        K1v = K1 + zeros(nVari-1,1);
        K2 = 0;%128*visco_slvw/pi/wellD^4*Qchara/Pcharaw*1;
        
        L = zeros(nVari-1,1);
        RHS(2*nAct+nQ) = sum(Q(nQ+1:nQ+nVari-1)) - Q0;
        Jac(2*nAct+nQ,2*nAct+nQ+1:2*nAct+nQ+nQ_const-1) = 1;
        
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
clear RHS
clear Jac
end


