function [Pres,DD] = Intial_Balance(nEle,AllEle,ConnList,P0,rate,DD0)
global Mat
P0 = P0 + Mat.Sxx*1e6;
Pbnd = P0(1);
H = Mat.h;
TipStates(1) = 1;
TipStates(nEle) = 1;
RHS = zeros(nEle*2,1);
EleL = AllEle(:,7);
Pres = P0;
DD = [zeros(nEle,1);-DD0];
Dn = DD(nEle+1:2*nEle);
%
eps = 1e-4;
dis = 1e6;
Jac = zeros(nEle*2,nEle*2);
while dis > eps
    % PP boundary condition -- PP = 0
    [~,CM] = BuildCoefMatix_Constant_pre(Mat,nEle,H,AllEle,TipStates);
    A110 = CM(1:nEle,1:nEle);
    A12 = CM(1:nEle,nEle+1:2*nEle);
    A21 = CM(nEle+1:2*nEle,1:nEle);
    A22 = CM(nEle+1:2*nEle,nEle+1:2*nEle);
    tempI = eye(nEle);
    BC = zeros(nEle*2,1);
    for ii = 1 : nEle
        BC(nEle+ii) =BC(nEle+ii) - Pres(ii);
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
    %Derivative of preCMure
    for ii = 1 : nEle
        Jac(ii,nEle+ii) = 1;
    end
    % Unit MPa, SI unit
    RHS(1:nEle) =A11*Dns-RDs;
    if dis == 1e6
        Dn = A11\RDs;
        RHS(1:nEle) =A11*Dn-RDs;
        DD(nEle+1:nEle*2) = Dn;
    end
    e = - DD(nEle+1:nEle*2);
    % Flow
    %Volumn
    V = e.*EleL*H;
    h = H;
    % A = e*h;
    disp('Building Jacob Matrix of Fluid Flow Equations');
    for i = 1: nEle
        nConn = ConnList(i,2);
        ki = e(i)^2/12;
        Li = AllEle(i,7);
        % Flow Flux Term
        for j = 1 : nConn
            if ConnList(i,j+2) > -0.1
                Conj =  ConnList(i,j+2);
                if (Conj) < 0.1
                    continue;
                end
                if e(Conj) < 1e-16
                    continue;
                end
                Lj = AllEle(Conj,7);
                A = (e(i)+e((Conj)))/2*h;
                dL = (Li + Lj)/2;
                kj =  e((Conj))^2/12;
                v = ki*Li+kj*Lj;
                kij = ki*kj*(Li+Lj)/(ki*Li+kj*Lj);
                u = ki*kj*(Li+Lj);
                
                dkj_dj = -e((Conj))/6;
                dA_dj = -h/2;
                dv_dj = dkj_dj*Lj;
                du_dj = dkj_dj*ki*(Li+Lj);
                T_difDj = dA_dj*kij + A*(du_dj*v - dv_dj*u)/v^2;
                T_difPj = 0;
                
                dki_di = -e(i)/6;
                dA_di = -h/2;
                du_di = dki_di*kj*(Li+Lj);
                dv_di = dki_di*Li;
                T_difDi = dA_di*kij + A*(du_di*v - dv_di*u)/v^2;
                T_difPi = 0;
                unitC = 1;
                Tij = A*kij/dL;
                %Upwind format
                if Pres(i) >= Pres((Conj))
                    b = Calc_Bw(Pres(i),1);
                    b_difp = Calc_Bw(Pres(i),2);
                    revmiu =Calc_miu(Pres(i),1);
                    revmiu_difp =Calc_miu(Pres(i),2);
                    Jac(nEle+i,nEle+i) = Jac(nEle+i,nEle+i) + Tij * ((b*revmiu_difp + revmiu*b_difp)*(Pres((Conj)) - Pres(i)) - revmiu*b)*unitC + T_difPi * (revmiu*b)*(Pres((Conj)) - Pres(i));
                    Jac(nEle+i,nEle+(Conj)) = Jac(nEle+i,nEle+(Conj)) +Tij * (revmiu*b)*unitC + T_difPj * (revmiu*b)*(Pres((Conj)) - Pres(i));
                else
                    b = Calc_Bw(Pres((Conj)),1);
                    b_difp = Calc_Bw(Pres((Conj)),2);
                    revmiu =Calc_miu(Pres((Conj)),1);
                    revmiu_difp =Calc_miu(Pres((Conj)),2);
                    Jac(nEle+i,nEle+i) = Jac(nEle+i,nEle+i) +  Tij * (-revmiu*b)*unitC + T_difPi * (revmiu*b)*(Pres((Conj)) - Pres(i));
                    Jac(nEle+i,nEle+(Conj)) = Jac(nEle+i,nEle+(Conj)) + Tij * (revmiu*b+(b*revmiu_difp + revmiu*b_difp)*(Pres((Conj)) - Pres(i)))*unitC + T_difPj * (revmiu*b)*(Pres((Conj)) - Pres(i));
                end
                if isnan(Jac(nEle+i,nEle+i))
                    keyboard;
                end
                Jac(nEle+i,i) = Jac(nEle+i,i) +T_difDi * (-1)* (revmiu*b)*(Pres((Conj)) - Pres(i));
                Jac(nEle+i,(Conj)) = Jac(nEle+i,(Conj)) +T_difDj * (-1)*(revmiu*b)*(Pres((Conj)) - Pres(i));
                RHS(i + nEle) = RHS(i + nEle) + Tij * (revmiu*b)*(Pres((Conj)) - Pres(i));
                if isnan(Jac(nEle+i,nEle+i))
                    keyboard;
                end
            end
        end
        % Inj
        if (i == nEle/2 )||( i == nEle/2+1)
            RHS(i+nEle) = RHS(i+nEle) +  rate;
        end
    end
    Jac(nEle+1,:) = 0;
    Jac(nEle*2,:) = 0;
    Jac(nEle+1,nEle+1) = 1;
    Jac(nEle*2,nEle*2) = 1;
    RHS(nEle+1) = 0;
    RHS(nEle*2) = 0;
    dX = -Jac\RHS;
    Pres = Pres+dX(nEle+1:nEle*2);
    DD(nEle+1:nEle*2) = DD(nEle+1:nEle*2) + dX(1:nEle);
    Dn = DD(nEle+1:2*nEle);
%
end
f = 1;
end