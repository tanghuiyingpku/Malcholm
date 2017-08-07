function  DD= CalcDD_Cross(dt,CurT,nEle,DD,Pres,AllEle,ConnList)
global Index IndexInv
global isMechActive;
global nwell well Mat Fluid;
global CpF_global XfF_global   InitialAperture
global Hbanking_global  ;
Schindex = zeros(nwell,1);
for welli = 1 : nwell
    for i = 1 : well{welli}.nSch
        if CurT >= well{welli}.Sch(i).t0 && CurT <= well{welli}.Sch(i).t1
            Schindex(welli) = i;
        end
    end
end
restart = 1;
e0 = -DD(nEle+1:2*nEle);
P_pre = Pres;
DD_pre = DD;

hbanking = zeros(nEle,Fluid.nprop);
Cp0 = zeros(nEle,Fluid.nprop);
Xf0 = zeros(nEle,Fluid.nfluid);
% Assume the maximum connection number is 4
Vel_sl = zeros(nEle,4);
visco_flv = zeros(nEle,1);
visco_slv = zeros(nEle,1);
dens_slv = zeros(nEle,1);
dens_slv0 = zeros(nEle,1);
dens_dpv = zeros(nEle,1);
densfv = zeros(nEle,Fluid.nfluid);
densf0v = zeros(nEle,Fluid.nfluid);
for i = 1 : nEle-1
    Cp0(i,:) = CpF_global(IndexInv(i),:);
    Xf0(i,:) = XfF_global(IndexInv(i),:);
    hbanking(i,:) = Hbanking_global(IndexInv(i),:);
end
Xf0(nEle,:) = Xf0(nEle-1,:)*0+1;
hbanking(nEle,:) = hbanking(nEle-1,:)*0;
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
%
InsituS = zeros(nEle*2,1);
for i = 1 : nEle
    %alpha = theta + 90;
    sinalp = AllEle(i,6);  %cosbet(i);
    cosalp = -AllEle(i,5); %-sinbet(i);
    InsituS(i) = (Sxx - Syy) * cosalp * sinalp - Sxy*(cosalp^2 - sinalp^2);
    InsituS(i + nEle) = Sxx * cosalp^2 + Syy * sinalp^2+2*Sxy*sinalp*cosalp;
end
BC_insitu = InsituS;
% Convergence Criteria
tol = 1e-5;
sigmaN = zeros(nEle,1);
ShearBC = zeros(nEle,1);
Cp = Cp0;
Xf = Xf0;
nQ = 0;
Q = zeros(50,1);
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
% Give initial opening to tip elements
DDcur = DD_pre;
for ii = 1 : nEle
    if abs(DDcur(ii+nEle) ) < InitialAperture
        DDcur(ii+nEle) = -InitialAperture;
    end
    %end
end
%
for ii = 1 : nEle
    isMechActive(ii) = 1;
end
isMechActive0 = isMechActive;
%Previous Time Step Values
iIterT = 0;
while restart > 0.1 || iIterT > 0.1
    %n Ds n Dn n Pressure
    isMechActive = isMechActive0;
    Jac0 = zeros(2*nEle,2*nEle);
    RHS = zeros(2*nEle,1);
    % RHS for solid equation set
    BC0   = zeros(2*nEle,1);
    % Project the in-situ stresses to fractures;
    restart = 0;
    DD = DDcur;
    P = P_pre/Pchara;
    % No fluid Element
    isRegroup = 1;
    % Regroup is for types changes
    Dn = DD(nEle+1:2*nEle);
    %Ds = DD(1:nEle);
    %while isRegroup > 0.1
    dis = 1e6;
    iIter = 0;
    BC = BC0;
    Jac = Jac0;
    while dis > tol || isRegroup > 0.1
        % Properties
        for ii = 1 : nEle
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
        if iIter < 1.1
            [extraBC,SS]  = BuildCoefMatix_Constant_test(nEle,Mat,AllEle);
        end
        disp('Coefficient Matrix Building finish');
        BC = BC + extraBC;
        % Elimination of Shear Displacements
        A110 = SS(1:nEle,1:nEle);
        A12 = SS(1:nEle,nEle+1:2*nEle);
        A21 = SS(nEle+1:2*nEle,1:nEle);
        A22 = SS(nEle+1:2*nEle,nEle+1:2*nEle);
        tempI = eye(nEle);
        
        EleL  = zeros(nEle,1);
        % For Pressure in Solid equation SF part Solid-Fluid
        e = -Dn;
        for ii = 1 : nEle
            EleL(ii) = AllEle(ii,7);
            if isMechActive(ii) < 0.1
                if isMechActive(ii) < -0.1
                    %  totally closed
                    BC(ii) = 0;
                    BC(nEle+ii) = 0;
                    A110(ii,:) = 0 ;
                    A110(ii,ii) = 1;
                    A12(ii,:) = 0 ;
                    A22(ii,:) = 0 ;
                    A22(ii,ii) = 1;
                    A21(ii,:) = 0 ;
                else
                    % Sliding, Only Ds is included in calculation
                    BC(nEle+ii) = 0;
                    A22(ii,:) = 0 ;
                    A22(ii,ii) = 1;
                    A21(ii,:) = 0 ;
                    BC(ii) = BC(ii)-BC_insitu(ii) + ShearBC(ii);
                end
            else
                BC(ii) = BC(ii)-BC_insitu(ii);
                BC(nEle+ii) =BC(nEle+ii) -BC_insitu(nEle+ii) - P(ii)*Pchara;
            end
        end
        invA11 = A110\tempI;
        R1 = BC(1:nEle);
        R2 = BC(nEle+1:2*nEle);
        A11 = -A21*invA11*A12 + A22;
        RD = R2 - A21*invA11*R1;
        %Only Include the Active Grids
        %Elements with fluid
        Jac(1:nEle,1:nEle) = A11;
        RDs = RD;
        Dns = Dn;
        Pact = P;
        %Derivative of pressure
        for ii = 1 : nEle
            if isMechActive(ii) > 0.1
                %-- Scaling--
                Jac(ii,nEle+ii) = 1*Pchara;
            end
        end
        % Unit MPa, SI unit
        RHS(1:nEle) =A11*Dns-RDs;
        %----------Scaling------------
        Jac(1:nEle,1:nEle) = Jac(1:nEle,1:nEle) * Dchara;
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
        for i = 1: nEle
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
                    Lj = AllEle(Index(Conj),7);
                    A = (e(i)+e(Index(Conj)))/2*h;
                    dL = (Li + Lj)/2;
                    kj = max(e(Index(Conj))^2/12,InitialAperture^2/12);
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
                        Jac(nEle+i,nEle+i) = Jac(nEle+i,nEle+i) +  Pchara*1/viscosl*(Tij * (dens_dp*(P(Index(Conj)) - P(i)) - dens_sl)*unitC + T_difPi * dens_sl*(P(Index(Conj)) - P(i)));
                        Jac(nEle+i,nEle+Index(Conj)) = Jac(nEle+i,nEle+Index(Conj)) + Pchara*1/viscosl * Tij * dens_sl*unitC + Pchara*T_difPj * dens_sl*(P(Index(Conj)) - P(i));
                        Vel_sl(i,j) = 1/viscosl * kij * (P(Index(Conj)) - P(i))/dL*Pchara;
                    else
                        viscosl =visco_slv(Index(Conj));
                        dens_sl = dens_slv(Index(Conj));
                        dens_dp = dens_dpv(Index(Conj));
                        Jac(nEle+i,nEle+i) = Jac(nEle+i,nEle+i) +  Pchara*1/viscosl * (Tij* (-dens_sl)*unitC + T_difPi * dens_sl*(P(Index(Conj)) - P(i)));
                        Jac(nEle+i,nEle+Index(Conj)) = Jac(nEle+i,nEle+Index(Conj)) +Pchara*1/viscosl*( Tij * (dens_sl+dens_dp*(P(Index(Conj)) - P(i)))*unitC + T_difPj * dens_sl*(P(Index(Conj)) - P(i)));
                        Vel_sl(i,j) = 1/viscosl * kij * (P(Index(Conj)) - P(i))/dL*Pchara;
                    end
                    if isnan(Jac(nEle+i,nEle+i))
                        keyboard;
                    end
                    Jac(nEle+i,i) = Jac(nEle+i,i) +1/viscosl*T_difDi * Dchara * (-1)* dens_sl*(P(Index(Conj)) - P(i))*Pchara;
                    Jac(nEle+i,Index(Conj)) = Jac(nEle+i,Index(Conj)) +1/viscosl*T_difDj* Dchara * (-1)*dens_sl*(P(Index(Conj)) - P(i))*Pchara;
                    RHS(i + nEle) = RHS(i + nEle) + 1/viscosl*Tij * dens_sl*(P(Index(Conj)) - P(i))*Pchara;
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
            Jac(nEle+i,nEle+i) = Jac(nEle+i,nEle+i) - (V(i)*dens_dp/dt*unitQ + (dwdp*dens_sl*h*EleL(i))/dt*unitQ)*1;
            if isnan(Jac(nEle+i,nEle+i))
                keyboard;
            end
            % For Opening Derivative
            Jac(nEle+i,i) = Jac(nEle+i,i) - dDdp*AllEle(i,7)*h*dens_sl/dt*unitQ*Dchara;
            %Newly Modified,Accumulation Terms
            dVdt = (V(i)*dens_sl-e0(i)*dens_sl0*h0*AllEle(i,7))/dt*unitQ;
            if abs(dVdt) < 1e-16
                dVdt = 0;
            end
            RHS(i+nEle) = RHS(i+nEle) - dVdt;
            %Carter's Leak
            %             RHS(i+nEle) = RHS(i+nEle) - CL*Li*h/sqrt(CurT+dt - Tau_global(IndexInv(i)))*2*dens_sl; % Two Faces Leak off
        end
        %% Well Part
        nQ = 0;
        for ii = 1 : nwell
            for jj = 1 :well{ii}.Sch(Schindex(ii)).nPf
                iele =  well{ii}.Sch(Schindex(ii)).Pf(jj);
                num = Index(well{ii}.Perfindex(iele,:));
                if min(num) > 1e-6
                    % Jacob for well
                    nQ = nQ + 1;
                    RHS(num(1)+nEle) = RHS(num(1)+nEle) + Q(nQ+1)/2*dens_sl(num(1));
                    RHS(num(2)+nEle) = RHS(num(2)+nEle) + Q(nQ+1)/2*dens_sl(num(1));
                end
            end
        end
        %%
        X0 = [Dns/Dchara;Pact];
        disp('Solve Linear Equation System');
        % Jac = DumpJacob(Jac);
        dx = -Jac\RHS;
        disp('Linear Equation Solved');
        X = X0 + dx;
        Dn = X(1:nEle)*Dchara;
        Pact  = X(nEle+1:nEle+nEle)*Pchara;
        P = Pact/Pchara;
        Ds = invA11*R1 - invA11*A12*Dn;
        for i = 1 : nEle
            if Dn(i) > -1e-7
                Dn(i) = 0;
            end
        end
        DD = [Ds;Dn];
        dis0 = dis;
        disp('   ');
        disD = norm(RHS(1:nEle));
        disP = norm(RHS(nEle+1:nEle+nEle))/1e1;
        xnormD = norm(dx(1:nEle))/(sum(abs(X(1:nEle))));
        xnormP = norm(dx(nEle:nEle+nEle))/(sum(abs(X(nEle:nEle+nEle))));
        % disStates = norm(isMechActive - isMechActive0);
        dis = max([disD,disP,xnormD,xnormP]);
        if abs(dis - dis0) < 1e-5
            dis = 0;
        end
        % dis = max([xnormD  xnormP xnormQ disQ]);
        fprintf('**Number of Iteration is : %d  Relative Convergence Norm Value times %f \n', iIter,dis);
        isRegroup = 0;
        disp('Check Elements Types...');
        disp('   ');
        fprintf('If Regroup ?: %d\n',isRegroup);
        %        if iIter > 5
        %       end
        %         if iIter >= 30
        %             dt = dt/3;
        %             restart =  1;
        %             break;
        %         end
    end
end
% check mass balance
%detQ = MassBalance(CurT,e0,e,Schindex,dens_slv0,dens_slv,h,dt,CL);
fprintf('**************\n');
Gtotals = 0;
Gtotaln = 0;

for i = 1 : nEle
    Dn = DD(i+nEle);
    Ds = DD(i);
    if DD(i+nEle) > 1e-8
        Dn = 0;
        Ds = 0 ;
    end
    Gtotals = Gtotals + 0.5*AllEle(i,7)*(Ds*(0-InsituS(i)) )*Mat.h;
    Gtotaln = Gtotaln + 0.5*AllEle(i,7)* (-P(i)-InsituS(i+nEle))*Dn*Mat.h;
end

dG = 0;
P = P * 1E6;
P_pre = P_pre*1e6;
for i = 1 : nEle
    Dn = DD(i+nEle);
    Ds = DD(i);
    if DD(i+nEle) > 1e-8
        Dn = 0;
        Ds = 0 ;
    end
    dG = dG + 0.5*AllEle(i,7)*((Ds-DD_pre(i))/dt*(0-InsituS(i)) + (Dn - DD_pre(i+nEle))/dt*(-P(i)-InsituS(i+nEle)))*Mat.h;
    if i < nEle
        dG = dG - 0.5*AllEle(i,7)*Dn *((-P(i)-InsituS(i+nEle))-(-P_pre(i)-InsituS(i+nEle)))/dt*Mat.h;
    else
        dG = dG - 0.5*AllEle(i,7)*Dn *((-P(i)-InsituS(i+nEle))-0)/dt*Mat.h;
    end
end

% 
% GI = 0;
% nQ = 0;
% for ii = 1 : nwell
%     wellOpen = 0;
%     for jj = 1 :well{ii}.Sch(Schindex(ii)).nPf
%         iele =  well{ii}.Sch(Schindex(ii)).Pf(jj);
%         num = Index(well{ii}.Perfindex(iele,:));
%         if min(num) > 1e-6
%             if wellOpen < 0.1
%                 nQ = nQ + 1;
%                 nQ = nQ + 1;
%                 GI = GI + Q(nQ)* (P(num(1))+InsituS(num(1)+nEle)) * dt;
%                 wellOpen = 1;
%             else
%                 nQ = nQ + 1;
%                 GI = GI + Q(nQ)* (P(num(1))+InsituS(num(1)+nEle)) * dt;
%             end
%         end
%     end
% end
% diss = 0;
% for i = 1 : nEle
%     isTip = 0;
%     for j = 1 : ConnList(i,2);
%         if ConnList(i,2+j) < 0.1
%             isTip = 1;
%         end
%     end
%     if isTip < 0.1 && ConnList(i,2) > 1.1
%         num1 =  ConnList(i,3) ;
%         num1 = Index(num1);
%         num2 =  ConnList(i,4);
%         num2 = Index(num2);
%         if num1 > 0.1 && num2 > 0.1
%             dP = (P(num1) - P(num2))/2/AllEle(i,7);
%             diss = diss + dP^2*AllEle(i,7)*abs(DD(i+nEle))^3/6/Mat.miu*Mat.h*dt;
%         end
%     else
%         num1 =  ConnList(i,3) ;
%         if num1 > 0.1
%             num1 = Index(num1);
%         end
%         num2 =  ConnList(i,4);
%         if num2 > 0.1
%             num2 = Index(num2);
%         end
%         if num1 > 0.1
%             dP = (P(num1) - P(i))/AllEle(i,7);
%             diss = diss + dP^2*AllEle(i,7)*abs(DD(i+nEle))^3/6/Mat.miu*Mat.h*dt;
%         end
%          if num2 > 0.1
%             dP = (P(num2) - P(i))/AllEle(i,7);
%             diss = diss + dP^2*AllEle(i,7)*abs(DD(i+nEle))^3/6/Mat.miu*Mat.h*dt;
%         end
%     end
% end
% fprintf('Max Pressure = %f\n',max(P));
% Gtotaln = GI- Gtotaln - diss;

end