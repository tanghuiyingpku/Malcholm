function dt =  FluidSolidFullCouple_Scaling_well(CurT,dt,isGrow)
% The boundary condition for fluid flow in no flow bounday
% In this function Shear Displacements has nothing to do with fluid flow so
% is not taken into direct calculation
% The elements with no fluid flow is not taken into calculation

%% Build Fluid Parts Jocob and RHS
global Index IndexInv nAct nAllEle_global
global isMechActive;
global TipStates EleType;
global nwell well Mat Fluid;
global  MaxEle ActEle;
global DD_global PresF_global AllEle_global ConnList_global;
global CpF_global XfF_global  Fractures;
global Hslurry_global Hbanking_global e_global Tau_global;
global InitialAperture;
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
Dsi = zeros(nAct,1);
DD_pre = zeros(nAct*2,1);
Cp0 = zeros(nAct,Fluid.nprop);
Xf0 = zeros(nAct,Fluid.nfluid);
Tipstate = zeros(nAct,1);
% Assume the maximum connection number is 4
Vel_sl = zeros(nAct,4);
Vel_fl = zeros(nAct*Fluid.nfluid,4);
Vel_prop = zeros(nAct*Fluid.nprop,4);
visco_flv = zeros(nAct,1);
visco_slv = zeros(nAct,1);
dens_slv = zeros(nAct,1);
dens_slv0 = zeros(nAct,1);
dens_dpv = zeros(nAct,1);
densfv = zeros(nAct,Fluid.nfluid);
densf0v = zeros(nAct,Fluid.nfluid);

for i = 1 : nAct
    e0(i) = -DD_global(MaxEle+IndexInv(i));
    Dsi(i) = -DD_global(IndexInv(i));
    DD_pre(i) = DD_global(IndexInv(i));
    DD_pre(i+nAct) = DD_global(MaxEle+IndexInv(i));
    P_pre(i) = PresF_global(IndexInv(i))*1e6;
    Cp0(i,:) = CpF_global(IndexInv(i),:);
    Xf0(i,:) = XfF_global(IndexInv(i),:);
    Tipstate(i) = TipStates(IndexInv(i),:);
    hslurry(i,:) = Hslurry_global(IndexInv(i),:);
    hbanking(i,:) = Hbanking_global(IndexInv(i),:);
end
hslurry0 = hslurry;
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
tol = 1e-5;
%How close the normal stress is calculated
epsd = 1e-8;
tolNorm = 1.05;
sigmaN = zeros(nAct,1);
ShearBC = zeros(nAct,1);
clear isMechActive0;
isMechActive0 = isMechActive;
Cp = Cp0;
Xf = Xf0;
nQ = 0;
Q = zeros(50,1);
%Size of well unkowns
for ii = 1 : nwell
    wellOpen = 0;
    for jj = 1 :well{ii}.Sch(Schindex(ii)).nPf
        iele =  well{ii}.Sch(Schindex(ii)).Pf(jj);
        num = Index(well{ii}.Perfindex(iele,:));
        if min(num) > 1e-6
            if wellOpen < 0.1
                nQ = nQ + 1;
                %% Jacob for well
                Q(nQ) = -Mat.Sxx*1.2*1e6/Pchara;
                wellOpen = 1;
                nQ = nQ + 1;
                Q(nQ) = well{ii}.Sch(Schindex(ii)).Pf_Q(jj,1)*2;
                %%
            else
                nQ = nQ + 1;
                Q(nQ) = well{ii}.Sch(Schindex(ii)).Pf_Q(jj,1)*2;
            end
            
        end
    end
end
Q = Q(1:nQ);
nQ_const = nQ;
for ii = 1 : nAct
    if abs(AllEle(ii,10) - 3)<1e-3
        % Closed Fracture
        xx= AllEle(ii,8)+epsd*AllEle(ii,5);
        yy = AllEle(ii,9)-epsd*AllEle(ii,6);
        [Sxxi,Syyi,Sxyi,~,~] = CalcPointStress_C(Mat,DD_pre,xx,yy,AllEle_global(1:nAllEle_global,:),ConnList_global(1:nAllEle_global,:));
        Sxxi = Sxxi + Mat.Sxx;
        Syyi = Syyi + Mat.Syy;
        Sxyi = Sxyi + Mat.Sxy;
        cosb = AllEle(ii,6);
        sinb = AllEle(ii,5);
        sigmaN(ii) = abs(Sxxi*sinb^2-2*Sxyi*sinb*cosb + Syyi*cosb^2)*1e6;
        if P_pre(ii) < sigmaN(ii)*tolNorm
            isMechActive(ii) = 0;
            DD_pre(ii+nAct) = 0;
            % Closed NF
            [e0(ii),~] = calcNfWidth(sigmaN(ii),P_pre(ii));
        else
            isMechActive(ii) = 1;
            e0(ii) = e0(ii) + calcNfWidth(0,0);
        end
    else
        isMechActive(ii) = 1;
    end
end
isMechActive0 = isMechActive;
%Previous Time Step Values
DDcur = DD_pre;
while restart > 0.1
    %n Ds n Dn n Pressure
    isMechActive = isMechActive0;
    Jac0 = zeros(2*nAct,2*nAct);
    RHS = zeros(2*nAct,1);
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
    isAdjustState= 0;
    while dis > tol || isRegroup > 0.1
        if norm(isMechActive-isMechActive0) > 0.1
            isAdjustState = 1;
        end
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
        % P < normal Stress not mechanically active or open
        e = -Dn;
        for ii = 1 : nAct
            if abs(AllEle(ii,10) - 3)<1e-3
                if isAdjustState > 0.1 && isMechActive(ii) < 0.1
                    xx= AllEle(ii,8)+epsd*AllEle(ii,5);
                    yy = AllEle(ii,9)-epsd*AllEle(ii,6);
                    [Sxxi,Syyi,Sxyi,~,~] = CalcPointStress_C(Mat,DD_pre,xx,yy,AllEle_global(1:nAllEle_global,:),ConnList_global(1:nAllEle_global,:));
                    Sxxi = Sxxi + Mat.Sxx;
                    Syyi = Syyi + Mat.Syy;
                    Sxyi = Sxyi + Mat.Sxy;
                    cosb = AllEle(ii,6);
                    sinb = AllEle(ii,5);
                    sigmaN(ii) = abs(Sxxi*sinb^2-2*Sxyi*sinb*cosb + Syyi*cosb^2)*1e6;
                end
                if isMechActive(ii) > 0.1
                    [temp,~] = calcNfWidth(0,0) ;
                    e(ii) = e(ii) + temp;
                else
                    [e(ii),~] = calcNfWidth( sigmaN(ii),P(ii)*Pchara) ;
                    
                end
            end
        end
        
        for ii = 1 : nAct
            EleL(ii) = AllEle(ii,7);
            if isMechActive(ii) < 0.1
                if isMechActive(ii) < -0.1
                    %  totally closed
                    BC(ii) = 0;
                    BC(nAct+ii) = 0;
                    A110(ii,:) = 0 ;
                    A110(ii,ii) = 1;
                    A12(ii,:) = 0 ;
                    A22(ii,:) = 0 ;
                    A22(ii,ii) = 1;
                    A21(ii,:) = 0 ;
                else
                    % Sliding, Only Ds is included in calculation
                    BC(nAct+ii) = 0;
                    A22(ii,:) = 0 ;
                    A22(ii,ii) = 1;
                    A21(ii,:) = 0 ;
                    BC(ii) = BC(ii)-BC_insitu(ii) + ShearBC(ii);
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
        for ii = 1 : nAct
            if isMechActive(ii) < 0.1 && isMechActive(ii) > -0.1
                RHS(ii) = 0;
            end
            if abs(AllEle_global(IndexInv(ii),10)-1) < 1e-3 && TipStates(IndexInv(ii)) < 0.1
                if e(ii) < InitialAperture*0.9
                    %Closed HF
                    RHS(ii) = 0;
                    Dns(ii) = -InitialAperture;
                    e(ii) = InitialAperture;
                    Jac(ii,:) = Jac(ii,:)*0;
                    Jac(ii,ii)=1;
                end
            end
        end
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
            h = Mat.h - sum(hbanking(i,:));
            h0 = Mat.h - sum(hbanking0(i,:));
            V(i) = e(i).*EleL(i)*h;
            %Flow
            nConn = ConnList(i,2);
            ki = e(i)^2/12;
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
                    kj =  e(Index(Conj))^2/12;
                    kij = (ki*Li+kj*Lj)/(Li+Lj);
                    if isMechActive(Index(Conj)) > 0.1
                        dkj_dj = -e(Index(Conj))/6;
                        dA_dj = -h/2;
                        T_difDj = dA_dj*kij + A*(dkj_dj*Lj)/(Li+Lj);
                        T_difPj = 0;
                    else
                        T_difDj =0;
                        [~,dwdp] = calcNfWidth( sigmaN(Index(Conj)),P(Index(Conj))*Pchara);
                        dwdp = dwdp * Pchara;
                        dkj_dPj = e(Index(Conj))/6*dwdp;
                        dA_dPj = dwdp*h/2;
                        T_difPj = dA_dPj*kij + A*(dkj_dPj*Lj)/(Li+Lj);
                    end
                    if isMechActive(i) > 0.1
                        dki_di = -e(i)/6;
                        dA_di = -h/2;
                        T_difDi = dA_di*kij + A*(dki_di*Li)/(Li+Lj);
                        T_difPi = 0;
                    else
                        T_difDi = 0;
                        [~,dwdp] = calcNfWidth(sigmaN(i),P(i)*Pchara);
                        dwdp = dwdp * Pchara;
                        dki_dPi = e(i)/6*dwdp;
                        dA_dPi = dwdp*h/2;
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
            
            if isMechActive(i) < 0.1
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
            %             CL =1e-6*5;
            %             if nAct > 20
            %                 RHS(i+nAct) = RHS(i+nAct) - dt*CL*Li*h/sqrt(CurT+dt - Tau_global(IndexInv(i)))*2*dens_sl; % Two Faces Leak off
            %             end
            %              Darcy Leakoff
            %Pore Pressure
            Pp = Mat.Pp;
            Leff = 100;
            % Matric Permeability /m^2
            krock =1e-15*1e-1;
            Jac(i+nAct,i+nAct) = Jac(i+nAct,i+nAct) - krock*Li*h/Leff;
            RHS(i+nAct) = RHS(i+nAct) - krock*Li*h/Leff*(P(i) - Pp);
        end
        %% Well Part
        nQ = 0;
        for ii = 1 : nwell
            for jj = 1 :well{ii}.Sch(Schindex(ii)).nPf
                iele =  well{ii}.Sch(Schindex(ii)).Pf(jj);
                num = Index(well{ii}.Perfindex(iele,:));
                dens_slurry = 1020;
               % dens_slurry = ConstructWell();
                RHS(num(1)+nAct) = RHS(num(1)+nAct) + well{ii}.Sch(Schindex(ii)).Pf_Q(jj,1)*dens_slurry;
                RHS(num(2)+nAct) = RHS(num(2)+nAct) + well{ii}.Sch(Schindex(ii)).Pf_Q(jj,1)*dens_slurry;
            end
        end
        
        
        X0 = [Dns/Dchara;Pact];
        disp('Solve Linear Equation System');
        % Jac = DumpJacob(Jac);
        dx = -Jac\RHS;
        disp('Linear Equation Solved');
        X = X0 + dx;
        Dn = X(1:nAct)*Dchara;
        Pact  = X(nAct+1:nAct+nAct)*Pchara;
        P = Pact/Pchara;
        Ds = invA11*R1 - invA11*A12*Dn;
        DD = [Ds;Dn];
        disp('   ');
        disD = norm(RHS(1:nAct));
        disP = norm(RHS(nAct+1:nAct+nAct))/1e1;
        xnormD = norm(dx(1:nAct))/(sum(abs(X(1:nAct))));
        xnormP = norm(dx(nAct:nAct+nAct))/(sum(abs(X(nAct:nAct+nAct))));
        % disStates = norm(isMechActive - isMechActive0);
        % dis = max([disD,disP,xnormD,xnormP,disStates]);
        dis = max([disD,disP xnormD  xnormP ]);
        fprintf('**Number of Iteration is : %d  Relative Convergence Norm Value times %f \n', iIter,dis);
        isRegroup = 0;
        disp('Check Elements Types...');
        disp('   ');
        fprintf('If Regroup ?: %d\n',isRegroup);
        %        if iIter > 5
        for im = 1 : nAct
            if Dn(im) > 1e-12 && isMechActive(im) > 0.1
                Dn(im) = 0;
                if abs(AllEle(im,10)-1)> 1.1
                    isMechActive(im) = 0;
                end
            end
        end
        %       end
        if iIter >= 40
            dt = dt/3;
            restart = 1;
            break;
        end
    end
    
end
%detQ = MassBalance(CurT,e0,e,Schindex,dens_slv0,dens_slv,h,dt,CL);
fprintf('**************\n');
%fprintf('Relative Mass Balance = : %f\n',detQ);
e = -Dn;
% %     for ii = 1 : nAct
% %         if abs(AllEle(ii,10) - 3)<1e-3
% %             if isMechActive(ii) > 0.1
% %                 [temp,~] = calcNfWidth(P(ii)) ;
% %                 e(ii) = e(ii) + temp;
% %             else
% %                 [e(ii),~] = calcNfWidth(P(ii)) ;
% %             end
% %         end
% %     end
% %Update Cp Need velocity
% % Calc Proppant Vel Vel_prop[nAct*nProp,4]
[Vprop_z,Vprop_x] = CalcProppantV(nAct,Cp,Xf,densfv,propdens,visco_slv,e,Vel_sl);
%
% %     % Calc Fluid Velocity Vel_fl
% %     Vfl = Vel_sl;% = CalcFluidV(nAct,dens_slv,Vel_sl,Vprop_x,Cp,propdens);
% %     %Update Xf
% %     nprop = Fluid.nprop;
Cp_new = UpdateCp_x(Cp0,nAct,Fluid.nprop,Mat.h,hbanking,hbanking0,ConnList,Vprop_x,Schindex,e,e0,dt,AllEle(:,7));
% %     eps2 = 0;%norm(Cp-Cp_new);
Cp = Cp_new;
% % %    Xf_new = UpdataXf(Mat.h,hbanking,hbanking0,densfv,densf0v,ConnList,Schindex,e,e0,nAct,Fluid.nfluid,Vfl,Xf0,Cp,dt,AllEle(:,7));
% %     Xf_new = Xf;
% %     eps1 = norm(Xf-Xf_new);
% %     Xf = Xf_new;
% %     if min(Xf(:,1)) < 0.99
% %         f = 1;
% %     end
% %     isOne = Xf_new(:,1)+ Xf_new(:,2);
% %     % Update Cp Due to Horizontal Moving, Settling acts like sink term dV =
% %     % Cp*Vz*dt*e*dx
% %     if Schindex(1) > 1.1
% %         f = 1;
% % %         ShowBanking(hbanking);
% % %         ShowCp(Cp);
% %     end
% %
% %     % Update Cp Due to settling
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
% end
disp('Current Time step Coupling Equation Converged');
% Give Value to global  Variables
for i = 1 : nAct
    PresF_global(IndexInv(i)) = P(i)*Pchara/1e6;
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
                
                wellOpen = 1;
            else
                nQ = nQ + 1;
                well{ii}.Sch(Schindex(ii)).Pf_Q(jj,1) = Q(nQ)/2;
                well{ii}.Sch(Schindex(ii)).Pf_Q(jj,2) = Q(nQ)/2;
                well{ii}.Sch(Schindex(ii)).PfQsl(jj,1) = Q(nQ)/2*dens_slurry;
                well{ii}.Sch(Schindex(ii)).PfQsl(jj,2) = Q(nQ)/2*dens_slurry;
            end
            
        end
    end
end
%%

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
        Pcharaw =   1;%abs(Mat.Sxx)*1e6;
        Kd = 0.65;
        Dperf = 18*1e-3*100/2.54;
        Dens = Dens*0.0083454045;
        n2 = 15;
        % unit dens: lbm/gal dp:in RATE bpm = 0.159/60 m^3/s
        Qunit = (1*60/0.159)^2;
        K1 = 0.2369*Dens/n2^2/Dperf^4/Kd^2*Qchara^2/Pcharaw*Qunit;
        % From psi to Pa
        if CurT > 1200
            K1 = K1 / 145 * 1E2 * 1E1 * 0.5;
        else
            K1 = K1 / 145 * 1E2 * 1E1;
        end
        K1v = K1 + zeros(nVari-1,1);
        K2 = 0;%128*visco_slvw/pi/wellD^4*Qchara/Pcharaw*1;
        
        L = zeros(nVari-1,1);
        RHS(2*nAct+nQ) = 0;%sum(Q(nQ+1:nQ+nVari-1)) - Q0;
        Jac(2*nAct+nQ,2*nAct+nQ) = 1;
        
        for ip = 1 : nVari-1
            RHS(2*nAct+nQ+ip) = 0;
            Jac(2*nAct+nQ+ip,2*nAct+nQ+ip) = 1;
        end
        %             if Q(nQ+ip) > 1e-12
        %                 K1i = K1v(ip);
        %                 num_perf = Index(well{ii}.Perfindex(Perf(ip),:));
        %                 RHS(2*nAct+nQ+ip) = Q(nQ)*Pchara -  P(num_perf(1))/Pcharaw*Pchara - K1i*Q(nQ+ip)^2/4;
        %                 %Derivative of Pressure in fracture
        %                 Jac(2*nAct+nQ +ip,nAct + num_perf(1)) =  -  1/2/Pcharaw*Pchara ;
        %                 Jac(2*nAct+nQ +ip,nAct + num_perf(2)) =  -  1/2/Pcharaw*Pchara ;
        %
        %                 %Derivative of Rate in fracture flow equation
        %                 Jac(nAct + num_perf(1),2*nAct+nQ +ip)   =  0.5 * dens_slwell;
        %                 Jac(nAct + num_perf(2),2*nAct+nQ +ip)   =  0.5 * dens_slwell;
        %
        %                 % Derivative of Rate in well equation
        %                 Jac(2*nAct+nQ +ip,2*nAct+nQ+ip) = Jac(2*nAct+nQ +ip,2*nAct+nQ +ip) - 2/4*K1i*Q(nQ+ip);
        %                 % Derivative of well head pressure in well equation
        %                 Jac(2*nAct+nQ +ip,2*nAct+nQ) = Pchara;
        %
        %                 if ip == 1
        %                     L(ip) = CalculateDis(well{ii}.heel', well{ii}.Perf(ip,:));
        %                     dPfric = K2*L(ip)*Q0/Qchara;
        %                 else
        %                     L(ip)  = CalculateDis(well{ii}.Perf(ip,:), well{ii}.Perf(ip-1,:));
        %                     dPfric = K2*L(1)*Q0/Qchara;
        %                     for jP = 1: ip-1
        %                         dPfric = dPfric + L(jP+1)*(Q0/Qchara - sum(Q(nQ+1:nQ+jP)))* K2;
        %                         Jac(2*nAct+nQ +ip, 2*nAct+nQ+jP) = Jac(2*nAct+nQ +ip, 2*nAct+nQ+jP)+ sum(L(jP+1:ip)) * K2;
        %                     end
        %                 end
        %                 RHS(2*nAct+nQ+ip) =  RHS(2*nAct+nQ+ip) - dPfric;
        %             else
        %                 %Q < 0
        %                 RHS(2*nAct+nQ+ip) = 0;
        %                 Jac(2*nAct+nQ+ip,2*nAct+nQ+ip) = 1;
        %             end
        %         end
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


