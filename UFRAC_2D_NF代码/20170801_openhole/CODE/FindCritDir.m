function [CritTheta,CritK] = FindCritDir(KI1,KI2)
global nTip TipStates Mat;
global IndexInv nAct 
global ConnList_global AllEle_global;
CritTheta = zeros(nTip,1);
CritK  = zeros(nTip,1);
%Find Growth Direction and Critical K calculation
for ii = 1 : nAct
    if TipStates(IndexInv(ii)) > 0.1 %&& TipStates(IndexInv(ii)) < 999
        %if ConnList(ii,3)*ConnList(ii,4) < -0.1 %&& KI1(TipStates(ii)) > 1e-6
            if abs(KI1(TipStates(IndexInv(ii)))) < 1e-8
                CritK(TipStates(IndexInv(ii))) = 0;
                CritTheta(TipStates(IndexInv(ii))) = 0;
                continue;
            end
            [type,~] = TipType(IndexInv(ii),AllEle_global,ConnList_global);
            a = -90;
            b = -30;
            [hasFind,c]= FindOpenTheta(a,b,KI1(TipStates(IndexInv(ii))),KI2(TipStates(IndexInv(ii))));
            while hasFind < 0.1
                a = a+60;
                b = b+60;
                [hasFind,c]= FindOpenTheta(a,b,KI1(TipStates(IndexInv(ii))),KI2(TipStates(IndexInv(ii))));
            end
            if hasFind < 0.1
                %没有找到临界角度
                error('Fail to Find the critical Angle');
            end
            CritK(TipStates(IndexInv(ii))) = KI1(TipStates(IndexInv(ii))) * 0.5*cosd(c/2)*(1+cosd(c))-KI2(TipStates(IndexInv(ii))) *3/2*cosd(c/2)*sind(c);
            if type > 1.9
                c = c + 180;
            end
            CritTheta(TipStates(IndexInv(ii))) = c;
        %else
         %   CritK(TipStates(ii)) = sign(KI1(TipStates(ii)))*sqrt(KI1(TipStates(ii))^2+KI2(TipStates(ii))^2);
       % end
    end
end
CritK =(CritK/Mat.K1HF).^2;
 
end
