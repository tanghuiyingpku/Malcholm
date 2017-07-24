function [CritTheta,CritK] = FindCritDir_Anisotropy2(KI1,KI2)
global nTip TipStates;
global IndexInv nAct
global ConnList_global AllEle_global;
global Mat;
CritTheta = zeros(nTip,1);
CritK  = zeros(nTip,1);
%Find Growth Direction and Critical K calculation
for ii = 1 : nAct
    if TipStates(IndexInv(ii)) > 0.1 && TipStates(IndexInv(ii)) < 998
        %if ConnList(ii,3)*ConnList(ii,4) < -0.1 %&& KI1(TipStates(ii)) > 1e-6
        if abs(KI1(TipStates(IndexInv(ii)))) < 1e-8
            CritK(TipStates(IndexInv(ii))) = 0;
            CritTheta(TipStates(IndexInv(ii))) = 0;
            continue;
        end
        [type,~] = TipType(IndexInv(ii),AllEle_global,ConnList_global);
        a = -90;
        b = 90;
        angle_ele = asind(AllEle_global(IndexInv(ii),5));
        [Kmax,angle]= FindMaxTheta(a,b,KI1(TipStates(IndexInv(ii))),KI2(TipStates(IndexInv(ii))),angle_ele,Mat.Kbeta,Mat.Kmax,Mat.Kmin);
        CritK(TipStates(IndexInv(ii))) = Kmax;
        if type > 1.9
            angle = angle + 180;
        end
        CritTheta(TipStates(IndexInv(ii))) = angle;
        %else
        %   CritK(TipStates(ii)) = sign(KI1(TipStates(ii)))*sqrt(KI1(TipStates(ii))^2+KI2(TipStates(ii))^2);
        % end
    end
end
    function [Kic,angle]= FindMaxTheta(a,b,KI1,KI2,angle_ele,beta,Kmax,Kmin)
        Kic = 0;
        angle = 0;
        for theta0 = a : 0.1: b
            k1 = 0.5*cosd(theta0/2)*(KI1*(1+cosd(theta0))-3*KI2*sind(theta0));
            k2 = 0.5*cosd(theta0/2)*(KI1*sind(theta0)+KI2*(3*cosd(theta0)-1));
            % Angle between Kmax and K
            alpha = beta - theta0 - angle_ele;
            if abs(alpha) < 1e-8
                Kalpha = Kmax;
            else
                kk = cotd(alpha);
                Kalpha = sqrt(1+kk^2)*sqrt(1/(kk^2/Kmax^2+1/Kmin^2));
            end
            Ki = (k1^2+k2^2)/Kalpha^2;
            if Ki > Kic
                Kic = Ki;
                angle = theta0;
            end
        end
    end
end
