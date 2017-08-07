function [e,sigmaN,ShearBC,isMechActive2] = updateStates(isMechActive,AllEle,P,nAct, epsd,Ds,Ds0,Dn,Mat, fric, Maxtau)
sigmaN = zeros(nAct,1);
e = zeros(nAct,1);
ShearBC = zeros(nAct,1);
for i = 1 : nAct
    % Judge fracture type
    if  AllEle(i,10) > 1.1
        [tau0,Sn0] =  CalcPointStress_C_BC_local_3(i,[Ds;Dn]);
        Sxxi = Mat.Sxx;
        Syyi = Mat.Syy;
        Sxyi = Mat.Sxy;
        cosb = AllEle(i,6);
        sinb = AllEle(i,5);
                InsituS = (Sxxi*sinb^2-2*Sxyi*sinb*cosb + Syyi*cosb^2+Sn0)*1e6;

        sigmaN(i) = -(Sxxi*sinb^2-2*Sxyi*sinb*cosb + Syyi*cosb^2+Sn0)*1e6;
        tau = abs((Sxxi- Syyi)*cosb*(-sinb)-Sxyi*(sinb^2 - cosb^2)+tau0)*1e6;
        if isMechActive(i) > 0.1
            [ei,~] = calcNfWidth_S(Ds0(i),0,0,InsituS);
            e(i) = -Dn(i) + ei;
           % continue;
        end
%         isMechActive(i) = 0;
%         if sigmaN(i) <  P(i)*1.1
%             isMechActive(i) = 1; 
%         end
%         if tau > fric*(sigmaN(i) - P(i))&&  isMechActive(i) < 0.1
%             isMechActive(i) =0;
%         end
%         if isMechActive(i) < -0.1
%             isMechActive(i) = 0;
%         end
%         if  tau <= fric*(sigmaN(i) - P(i))*1.1 && isMechActive(i) < 0.1
%             isMechActive(i) = -2;
%         end
        % allele(i,10) = 5: First opened NF
        if isMechActive(i) < 0.1 % %sliding element
            ShearBC(i) = abs(fric*(sigmaN(i) - P(i)));
            if ShearBC(i) > Maxtau
                ShearBC(i) = Maxtau;
            end
            if sign(Ds(i)-Ds0(i)) > 0.1;
                ShearBC(i) = -ShearBC(i);
            end
        end
        if isMechActive(i) < 1e-5
            [ei,~] = calcNfWidth_S( Ds0(i),sigmaN(i),P(i),InsituS);
            e(i) = ei;
        else
            [ei,~] = calcNfWidth_S(Ds0(i),0,0,InsituS);
            e(i) = -Dn(i) + ei;
        end
    else
        isMechActive(i) = 1;
        e(i)  = -Dn(i);
        if e(i) <= 1e-20
            e(i) = 1e-16;
        end
    end
end
isMechActive2 = isMechActive;

end