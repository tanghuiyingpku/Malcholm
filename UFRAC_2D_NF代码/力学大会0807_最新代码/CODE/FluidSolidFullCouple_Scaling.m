function dt =  FluidSolidFullCouple_Scaling(CurT,dt,isGrow)
% The boundary condition for fluid flow in no flow bounday
% In this function Shear Displacements has nothing to do with fluid flow so
% is not taken into direct calculation
% The elements with no fluid flow is not taken into calculation
global Index IndexInv nAct nAllEle_global
global isMechActive;
global TipStates EleType;
global nwell well Mat Fluid;
global InitialAperture MaxEle;
global DD_global PresF_global AllEle_global ConnList_global;
global CpF_global XfF_global  Fractures;
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
for ii = nAct:-1:1
    if abs(Tipstate(ii)) > 1e-8
        labor = ConnList(ii,3);
        if labor < 0.1
            labor = ConnList(ii,4);
        end
        labor = Index(labor);
        if abs((DDcur(ii+nAct))+InitialAperture) < 1e-8
%             DDcur(ii+nAct) = DD_pre(labor+nAct)/sqrt(2);
%             P_pre(ii) = P_pre(labor)/2;
        end
    end
end
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
    %   Ds = DD(1:nAct);
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
        %dens_slv = dens_slv0;
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
        end
        h = Mat.h;
        %Pore Pressure
        Pp = Mat.Pp*1e6/Pchara;
        % Effective LeakOff Length\
        Leff = 100;
        % Matric Permeability /m^2
        krock = 0;%1e-15*1e-1;
        %Normal Opening :Open--> is negative
        if min(e) < 1e-16
            %  keyboard;
            ('Negative Opening');
        end
        %Volumn
        V = e.*EleL*Mat.h;
        % A = e*h;
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
            CL = 5*1e-7;
            if nAct > 41
                RHS(i+nAct) = RHS(i+nAct) - CL*Li*h/sqrt(CurT+dt - Tau_global(IndexInv(i)))*2*dens_sl; % Two Faces Leak off 
            end
%             if TipStates(IndexInv(i)) > 0.1
%                 Jac(i+nAct,:) = 0;
%                 Jac(i+nAct,i+nAct) = 1;
%                 RHS(i+nAct) = 0;
%                 Pact(i) = -Mat.Sxx;
%             end
        end
        for ii = 1 : nwell
            for jj = 1 :well{ii}.Sch(Schindex(ii)).nPf
                iele =  well{ii}.Sch(Schindex(ii)).Pf(jj);
                num = Index(well{ii}.Perfindex(iele,:));
                if min(num) > 1e-6
                    RHS(num(1)+nAct) = RHS(num(1)+nAct) + well{ii}.Sch(Schindex(ii)).PfQsl(jj,1);
                    RHS(num(2)+nAct) = RHS(num(2)+nAct) + well{ii}.Sch(Schindex(ii)).PfQsl(jj,2);
                end
            end
        end
        
        X0 = [Dns/Dchara;Pact];
        disp('Solve Linear Equation System');
        Jac(1:nAct,:) = Jac(1:nAct,:);
        RHS(1:nAct) = RHS(1:nAct);
        Jac = DumpJacob(Jac);
        dx = -Jac\RHS;
        disp('Linear Equation Solved');
        X = X0 + dx;
        Dn = X(1:nAct)*Dchara;
        
        if restart == 1
            break;
        end
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
        dis = max([disD,disP xnormD  xnormP]);
        fprintf('**Number of Iteration is : %d  Relative Convergence Norm Value times %f \n', iIter,dis);
        isRegroup = 0;
        disp('Check Elements Types...');
        disp('   ');
        fprintf('If Regroup ?: %d\n',isRegroup);
%        if iIter > 5
            for im = 1 : nAct
                if Dn(im) > -1e-12 && isMechActive(im) > 0.1
                    Dn(im) = 0;
                    if abs(AllEle(im,10)-1)> 1.1
                        isMechActive(im) = 0;
                    end
                end
            end
 %       end
        if iIter >= 15
            dt = dt/3;
            restart = 1;
            break;
        end
    end
end
detQ = MassBalance(CurT,e0,e,Schindex,dens_slv0,dens_slv,h,dt,CL);
fprintf('**************\n');
fprintf('Relative Mass Balance = : %f\n',detQ);

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
end


