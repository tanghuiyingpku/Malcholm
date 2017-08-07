function [isMechActive,ShearBC,sigmaN,e0,tau] = CheckStates_backup(nAct, AllEle, epsd,DD_pre,P_pre,Mat, fric, Maxtau)
global isMechActive_global;

sigmaN = zeros(nAct,1);
ShearBC = zeros(nAct,1);
tau = zeros(nAct,1);
isMechActive = isMechActive_global(1:nAct);
e0 = zeros(nAct,1);
for i = 1 : nAct
    % Judge fracture type
    if  AllEle(i,10) > 1.1
        [tau0,Sn0] =  CalcPointStress_C_BC_local_3(i,DD_pre);
        Sxxi = Mat.Sxx;
        Syyi = Mat.Syy;
        Sxyi = Mat.Sxy;
        cosb = AllEle(i,6);
        sinb = AllEle(i,5);
        InsituS = (Sxxi*sinb^2-2*Sxyi*sinb*cosb + Syyi*cosb^2+Sn0)*1e6;
        sigmaN(i) = -(Sxxi*sinb^2-2*Sxyi*sinb*cosb + Syyi*cosb^2+Sn0)*1e6;
        tau(i) = ((Sxxi- Syyi)*cosb*(-sinb)-Sxyi*(sinb^2 - cosb^2)+tau0)*1e6;
        if isMechActive(i) < 1e-5
            if sigmaN(i) < P_pre(i)*1.1
                isMechActive(i) = 1;
            end
            if sigmaN(i) >= P_pre(i)*1.1
                isMechActive(i) =0;
            end
            if  abs(tau(i)) <= fric*(sigmaN(i) - P_pre(i))*1.1 && isMechActive(i) < 0.1
                isMechActive(i) = -2;
            end
            if isMechActive(i) < 0.1 %sliding element
                %isMechActive(i) = 0;
                % Shear stress larger than maximum
                ShearBC(i) = abs(fric*(sigmaN(i) - P_pre(i)));
                if ShearBC(i) > Maxtau
                    ShearBC(i) = Maxtau;
                end
                if sign(tau(i)) < 0.1
                    ShearBC(i) = -ShearBC(i)*1;
                end
                if abs(DD_pre(i)) < 0.1
                    ShearBC(i) = 0;
                end
                % Default equals to 0
                %  ShearBC(i) = 0;
            end
            
            %closed and sliding element
            [ei,~] = calcNfWidth_S( DD_pre(i),sigmaN(i),P_pre(i),InsituS);
            e0(i) = ei;
        else
            %open NF element
            [ei,~] = calcNfWidth_S(DD_pre(i),0,0,InsituS);
            e0(i) = -DD_pre(i+nAct) + ei;
        end
    else
        isMechActive(i) = 1;
        e0(i)  = -DD_pre(i+nAct);
    end
end
end