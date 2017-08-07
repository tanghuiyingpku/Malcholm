function  [dt,DD_new,e_new,Pre_new,AllElenew] =FluidSolidFullCouple_Noflow_Scaling(CurT,EleType,Mat,P_pre,DD_pre,dt,nAllEle,AllEle,ConnList,nFracture,Fractures,isGrow)
% The boundary condition for fluid flow in no flow bounday
% In this function Shear Displacements has nothing to do with fluid flow so
% is not taken into direct calculation
% The elements with no fluid flow is not taken into calculation
global Index IndexInv nAct
global isMechActive;
global TipStates;
global nwell well;
global InitialAperture;
clear isMechActive0;
restart = 1;
e0 = DD_pre(nAct+1:2*nAct)*(-1);
P_pre = P_pre*1.02*1e6;
unitQ = 1;
unitC = 1;
Sxx = Mat.Sxx*1e6;
Syy = Mat.Syy*1e6;
Sxy = Mat.Sxy*1e6;
%Scaling characteristic value
Pchara = Sxx;
Dchara = 1e-3;
Vchara = 1e-3*Mat.h*0.1;
%
InsituS = zeros(nAct*2,1);
for i = 1 : nAct
    %alpha = theta + 90;
    sinalp = AllEle(IndexInv(i),6);  %cosbet(i);
    cosalp = -AllEle(IndexInv(i),5); %-sinbet(i);
    InsituS(i) = (Sxx - Syy) * cosalp * sinalp - Sxy*(cosalp^2 - sinalp^2);
    InsituS(i + nAct) = Sxx * cosalp^2 + Syy * sinalp^2+2*Sxy*sinalp*cosalp;
end
BC_insitu = InsituS;
% Convergence Criteria
tol = 1e-4;
%How close the normal stress is calculated
epsd = 1e-8;
tolNorm = 1.2;
Dsi = DD_pre(1:nAct);
sigmaN = zeros(nAct,1);
ShearBC = zeros(nAct,1);
clear isMechActive0;
isMechActive0 = isMechActive;

for ii = 1 : nAct
    if abs(AllEle(IndexInv(ii),10) - 3)<1e-3
        if isMechActive(ii) > 0.1
            e0(ii) = e0(ii) + calcNfWidth(P_pre(ii));
        else
            % Closed Fracture
            xx= AllEle(IndexInv(ii),8)+epsd*AllEle(IndexInv(ii),5);
            yy = AllEle(IndexInv(ii),9)-epsd*AllEle(IndexInv(ii),6);
            [Sxxi,Syyi,Sxyi,~,~] = CalcPointStress_C(Mat,DD_pre,xx,yy,AllEle,ConnList);
            Sxxi = Sxxi + Mat.Sxx;
            Syyi = Syyi + Mat.Syy;
            Sxyi = Sxyi + Mat.Sxy;
            cosb = AllEle(IndexInv(ii),6);
            sinb = AllEle(IndexInv(ii),5);
            sigmaN(ii) = abs(Sxxi*sinb^2-2*Sxyi*sinb*cosb + Syyi*cosb^2)*1e6;
            if P_pre(ii) < sigmaN(ii)*tolNorm
                DD_pre(ii+nAct) = 0;
                % Closed NF
                [e0(ii),~] = calcNfWidth(P_pre(ii));
            else
                isMechActive(ii) = 1;
                e0(ii) = e0(ii) + calcNfWidth(P_pre(ii));
            end
        end
    else
        isMechActive(ii) = 1;
    end
end
aa = size(isMechActive0);
aa2 = size(isMechActive);
if aa~=aa2
    f = 1;
end
if sum(isMechActive0 - isMechActive) > 0.1
    dif = 1;
end
%Previous Time Step Values
DDcur = DD_pre;
for ii = nAct:-1:1
    if abs(TipStates(IndexInv(ii))) > 1e-8
        labor = ConnList(IndexInv(ii),3);
        if labor < 0.1
            labor = ConnList(IndexInv(ii),4);
        end
        labor = Index(labor);
        if abs((DDcur(ii+nAct))+InitialAperture) < 1e-8
            DDcur(ii+nAct) = DD_pre(labor+nAct)/sqrt(3);
        end
    end
end
while restart > 0.1
    %n Ds n Dn n Pressure
    Jac0 = zeros(2*nAct,2*nAct);
    RHS = zeros(2*nAct,1);
    % RHS for solid equation set
    BC0   = zeros(2*nAct,1);
    % Project the in-situ stresses to fractures;
    restart = 0;
    DD = DDcur;
    P = P_pre;
    % No fluid Element
    isRegroup = 1;
    % Regroup is for types changes
    Dn = DD(nAct+1:2*nAct);
    Ds = DD(1:nAct);
    %while isRegroup > 0.1
    dis = 1e6;
    iIter = 0;
    BC = BC0;
    Jac = Jac0;
    while dis > tol || isRegroup > 0.1
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
                [extraBC,SS]  = BuildCoefMatix_Constant(Mat,nFracture,Fractures,nAllEle,AllEle,ConnList,1);
            else
                [extraBC,SS]  =BuildCoefMatix_Constant(Mat,nFracture,Fractures,nAllEle,AllEle,ConnList,0);
            end
        else
            if iIter < 1.1 && isGrow > 0.1
                [extraBC,SS]  = BuildCoefMatix_H2(Mat,nFracture,Fractures,nAllEle,AllEle,ConnList,1);
            else
                if iIter < 1.1
                    [extraBC,SS]  = BuildCoefMatix_H2(Mat,nFracture,Fractures,nAllEle,AllEle,ConnList,0);
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
            if abs(AllEle(IndexInv(ii),10) - 3)<1e-3
                if isMechActive(ii) > 0.1
                    [temp,~] = calcNfWidth(P(ii)) ;
                    e(ii) = e(ii) + temp;
                else
                    [e(ii),~] = calcNfWidth(P(ii)) ;
                end
            end
        end
        
        
        for ii = 1 : nAct
            EleL(ii) = AllEle(IndexInv(ii),7);
            if isMechActive(ii) < 0.1
                if isMechActive(ii) > -0.1
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
                BC(nAct+ii) =BC(nAct+ii) -BC_insitu(nAct+ii) - P(ii);
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
        Pp = Mat.Pp;
        % Effective LeakOff Length\
        Leff = 100;
        % Matric Permeability /m^2
        krock =0;% 1e-15*1e-3;
        %Normal Opening :Open--> is negative
        if min(e) < 1e-16
            %  keyboard;
            ('Negative Opening');
        end
        %Volumn
        V = e.*EleL*h;
        
        % A = e*h;
        disp('Building Jacob Matrix of Fluid Flow Equations');
        for i = 1: nAct
            %Flow
            nConn = ConnList(IndexInv(i),2);
            ki = e(i)^2/12;
            Li = AllEle(IndexInv(i),7);
            % Flow Flux Term
            for j = 1 : nConn
                if ConnList(IndexInv(i),j+2) > -0.1
                    Conj =  ConnList(IndexInv(i),j+2);
                    if Index(Conj) < 0.1
                        continue;
                    end
                    Lj = AllEle(Conj,7);
                    A = (e(i)+e(Index(Conj)))/2*h;
                    dL = (Li + Lj)/2;
                    kj =  e(Index(Conj))^2/12;
                    v = ki*Li+kj*Lj;
                    kij = ki*kj*(Li+Lj)/(ki*Li+kj*Lj);
                    u = ki*kj*(Li+Lj);
                    
                    if isMechActive(Index(Conj)) > 0.1
                        dkj_dj = -e(Index(Conj))/6;
                        dA_dj = -h/2;
                        dv_dj = dkj_dj*Lj;
                        du_dj = dkj_dj*ki*(Li+Lj);
                        T_difDj = dA_dj*kij + A*(du_dj*v - dv_dj*u)/v^2;
                        T_difPj = 0;
                    else
                        T_difDj =0;
                        [~,dwdp] = calcNfWidth(P(Index(Conj)));
                        dkj_dPj = e(Index(Conj))/6*dwdp;
                        dA_dPj = dwdp*h/2;
                        dv_dPj = dkj_dPj*Lj;
                        du_dPj = dkj_dPj*ki*(Li+Lj);
                        T_difPj = dA_dPj*kij + A*(du_dPj*v - dv_dPj*u)/v^2;
                    end
                    if isMechActive(i) > 0.1
                        dki_di = -e(i)/6;
                        dA_di = -h/2;
                        du_di = dki_di*kj*(Li+Lj);
                        dv_di = dki_di*Li;
                        T_difDi = dA_di*kij + A*(du_di*v - dv_di*u)/v^2;
                        T_difPi = 0;
                    else
                        T_difDi = 0;
                        [~,dwdp] = calcNfWidth(P(i));
                        dki_dPi = e(i)/6*dwdp;
                        dA_dPi = dwdp*h/2;
                        dv_dPi = dki_dPi*Li;
                        du_dPi = dki_dPi*kj*(Li+Lj);
                        T_difPi = dA_dPi*kij + A*(du_dPi*v - dv_dPi*u)/v^2;
                    end
                    Tij = A*kij/dL;
                    %Upwind format
                    if P(i) >= P(Index(Conj))
                        b = Calc_Bw(P(i),1);
                        b_difp = Calc_Bw(P(i),2);
                        revmiu =Calc_miu(P(i),1);
                        revmiu_difp =Calc_miu(P(i),2);
                        Jac(nAct+i,nAct+i) = Jac(nAct+i,nAct+i) + Pchara * Tij * ((b*revmiu_difp + revmiu*b_difp)*(P(Index(Conj)) - P(i)) - revmiu*b)*unitC + T_difPi * (revmiu*b)*(P(Index(Conj)) - P(i));
                        Jac(nAct+i,nAct+Index(Conj)) = Jac(nAct+i,nAct+Index(Conj)) + Pchara * Tij * (revmiu*b)*unitC + T_difPj * (revmiu*b)*(P(Index(Conj)) - P(i));
                    else
                        b = Calc_Bw(P(Index(Conj)),1);
                        b_difp = Calc_Bw(P(Index(Conj)),2);
                        revmiu =Calc_miu(P(Index(Conj)),1);
                        revmiu_difp =Calc_miu(P(Index(Conj)),2);
                        Jac(nAct+i,nAct+i) = Jac(nAct+i,nAct+i) +  Tij*Pchara * (-revmiu*b)*unitC + T_difPi * (revmiu*b)*(P(Index(Conj)) - P(i));
                        Jac(nAct+i,nAct+Index(Conj)) = Jac(nAct+i,nAct+Index(Conj)) + Tij *Pchara* (revmiu*b+(b*revmiu_difp + revmiu*b_difp)*(P(Index(Conj)) - P(i)))*unitC + T_difPj * (revmiu*b)*(P(Index(Conj)) - P(i));
                    end
                    if isnan(Jac(nAct+i,nAct+i))
                        keyboard;
                    end
                    Jac(nAct+i,i) = Jac(nAct+i,i) +T_difDi * Dchara * (-1)* (revmiu*b)*(P(Index(Conj)) - P(i));
                    Jac(nAct+i,Index(Conj)) = Jac(nAct+i,Index(Conj)) +T_difDj* Dchara * (-1)*(revmiu*b)*(P(Index(Conj)) - P(i));
                    RHS(i + nAct) = RHS(i + nAct) + Tij * (revmiu*b)*(P(Index(Conj)) - P(i));
                    if isnan(Jac(nAct+i,nAct+i))
                        keyboard;
                    end
                end
            end
            b = Calc_Bw(P(i),1);
            b_difp = Calc_Bw(P(i),2);
            b_init = Calc_Bw(P_pre(i),1);
            if isMechActive(i) < 0.1
                [~,dwdp] = calcNfWidth(P(i));
                dDdp = 0;
            else
                dwdp = 0;
                dDdp = -1;
            end
            % Accum Inj and Leak
            % Accum
            % For Pressure Derivative
            Jac(nAct+i,nAct+i) = Jac(nAct+i,nAct+i) - (V(i)*b_difp/dt*unitQ + (dwdp*b*h*EleL(i))/dt*unitQ)*Pchara;
            if isnan(Jac(nAct+i,nAct+i))
                keyboard;
            end
            % For Opening Derivative
            Jac(nAct+i,i) = Jac(nAct+i,i) - dDdp*AllEle(IndexInv(i),7)*h*b/dt*unitQ*Dchara;
            %Newly Modified,Accumulation Terms
            dVdt = (V(i)*b-e0(i)*b_init*h*AllEle(IndexInv(i),7))/dt*unitQ;
            if abs(dVdt) < 1e-16
                dVdt = 0;
            end
            RHS(i+nAct) = RHS(i+nAct) - dVdt;
            %Leak
            Jac(i+nAct,i+nAct) = Jac(i+nAct,i+nAct) - krock*Li*h/Leff;
            RHS(i+nAct) = RHS(i+nAct) - krock*Li*h/Leff*(P(i) - Pp);
        end
        
        % Horizontal Well
       HorizontalWellModule(CurT,P)
       
        for ii = 1 : nwell
            for jj = 1 : well{ii}.nPerf
                num = IndexInv(well{ii}.Perfindex(jj));
                RHS(num+nAct) = RHS(num+nAct) + well{ii}.PerfQ(jj);
            end
        end
        
        X0 = [Dns/Dchara;Pact/Pchara];
        disp('Solve Linear Equation System');
        Jac(1:nAct,:) = Jac(1:nAct,:) / Pchara;
        Jac(nAct+1:2*nAct,:) = Jac(nAct+1:2*nAct,:) / Vchara;
        RHS(1:nAct) = RHS(1:nAct) / Pchara;
        RHS(nAct+1:2*nAct) = RHS(nAct+1:2*nAct)/Vchara;
        dx = -Jac\RHS;
        
        disp('Linear Equation Solved');
        X = X0 + dx;
        Dn = X(1:nAct)*Dchara;
        Pact  = X(nAct+1:nAct+nAct)*Pchara;
        P = Pact;
        if isnan(P(1))
            keyboard;
        end
        Ds = invA11*R1 - invA11*A12*Dn;
        disp('   ');
        disD = norm(RHS(1:nAct));
        disP = norm(RHS(nAct+1:nAct+nAct));
        xnormD = norm(dx(1:nAct))/(sum(abs(X(1:nAct))));
        xnormP = norm(dx(nAct:nAct+nAct))/(sum(abs(X(nAct:nAct+nAct))));
        disStates = norm(isMechActive - isMechActive0);
        dis = max([disD,disP,xnormD,xnormP,disStates]);
        fprintf('**Number of Iteration is : %d  Relative Convergence Norm Value times %f \n', iIter,dis);
        DD = [Ds;Dn];
        
        isRegroup = 0;
        disp('Check Elements Types...');
        disp('   ');
        fprintf('If Regroup ?: %d\n',isRegroup);
        if iIter >= 30
            dt = dt/2;
            restart = 1;
            break;
        end
    end

end
disp('Current Time step Coupling Equation Converged');
DD_new = DD;
e_new = e;
Pre_new = P/1e6;
AllElenew = AllEle;
end


