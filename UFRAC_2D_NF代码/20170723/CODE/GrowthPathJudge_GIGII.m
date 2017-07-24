function  GrowthPathJudge_GIGII(nf0,CritTheta,K1c,K2c)
disp('** Judge the intersection behavior with Virtual Element Method')
% global variables
global  Tipcoordinate TipStatesInv  Mat  XfF_global Fractures
global TipStates ;
global nActTip nTip;
global IndexInv nAct Index MaxEle
global DD_global PresF_global AllEle_global  ConnList_global nAllEle_global
global HasInter

for i = 1 : nActTip
    nf = TipStatesInv(i);
    if TipStates(nf) > 998 && nf == nf0
        K1 = K1c(i);
        K2 = K2c(i);
        [Tiptype,~] = TipType(nf,AllEle_global,ConnList_global);
        type0 = ConnList_global(nf,1);
        con = ConnList_global(nf,3:5);
        nfcon = zeros(2,1);
        count = 1;
        cpoint = [AllEle_global(nf,3),AllEle_global(nf,4)];
        angleH = asind(AllEle_global(nf,5));
        if Tiptype > 1.1
            angleH = -asind(AllEle_global(nf,5));
            cpoint = [AllEle_global(nf,1),AllEle_global(nf,2)];
        end
        for j = 1 : 3
            typej = ConnList_global(con(j),1);
            if typej ~= type0
                nfcon(count) = con(j);%global numbering
                count = count + 1;
            end
        end
        Pinit = PresF_global(nf);
        Xfinit = XfF_global(nf,:);
        initAngle = CritTheta(i);
        
        %% Crossing Criteria
        rc = zeros(2,1);
        ds = zeros(2,1);
        for j = 1 : 2
            nb = nfcon(j);
            sinnf = AllEle_global(nb,5);
            alpha = asind(sinnf);
            vecnf = [AllEle_global(nb,8)-cpoint(1),AllEle_global(nb,9)-cpoint(2)];
            %[Cos Sin]
            veccoord = [AllEle_global(nb,6) AllEle_global(nb,5)];
            if Tiptype < 1.1
                theta0 = alpha - angleH;
            else
                theta0 = alpha - angleH + 180;
            end
            % Whether along the local coordiate
            vecdir = dot(vecnf,veccoord);
            if vecdir < -1e-10
                theta0 = theta0 + 180;
            end
            id = ConnList_global(nfcon(1));
            Fang = Fractures{id}.Fang;
            So = Fractures{id}.Coh;
            To = 2;
            if abs(theta0) > 180
                if theta0 > 0.1
                    theta0 = theta0 - 360;
                else
                    theta0 = theta0 + 360;
                end
            end
            [~,rc(j)] = Criteria(angleH,Mat.Sxx,Mat.Syy,Mat.Sxy,theta0,K1,K2,To);
            ds(j) = M_C_stability(angleH,Mat.Sxx,Mat.Syy,Mat.Sxy,theta0,K1,K2,rc(j),Fang,So);
        end
        
        if min(ds) < -1e-16
            % Sliding
            isCross = 0 ;
        else
            isCross = 1;
        end
        %%
        disp('Dilate into NF 1');
        
        if abs(DD_global(MaxEle+nfcon(1))) < 1e-7 %%&& abs(DD_global(nfcon(1)))< 1e-7 
            PresF_global(nfcon(1)) = Pinit;
            DD_global(MaxEle+nfcon(1)) = 0;
            XfF_global(nfcon(1),:) = Xfinit;
            AllEle_global(nfcon(1),10) = 3;
            point1 = AllEle_global(nfcon(1),1:2);
            dis =  sqrt((point1(1) - cpoint(1))^2 + (point1(2) - cpoint(2))^2 );
            if dis < 1e-6
                Tipcoordinate(nfcon(1),1:2) = AllEle_global(nfcon(1),3:4);
            else
                Tipcoordinate(nfcon(1),1:2) = point1;
            end
            TipStates(nfcon(1)) = nTip+1;
            Index(nfcon(1)) = nAct+1;
            IndexInv(nAct+1) = nfcon(1);
            nAct = nAct + 1;
        else
            if Index(nfcon(1)) < 0.1
                Index(nfcon(1)) = nAct+1;
                IndexInv(nAct+1) = nfcon(1);
                nAct = nAct + 1;
            end
        end
        if  abs(DD_global(MaxEle+nfcon(2))) < 1e-7%% && abs(DD_global(nfcon(2)))< 1e-7 
            PresF_global(nfcon(2)) = Pinit;
            DD_global(MaxEle+nfcon(2)) = 0;
            nTip = nTip + 1;
            disp('Dilate into NF 2');
            AllEle_global(nfcon(2),10) = 3;
            point1 = AllEle_global(nfcon(2),1:2);
            dis =  sqrt((point1(1) - cpoint(1))^2 + (point1(2) - cpoint(2))^2 );
            if dis < 1e-6
                Tipcoordinate(nfcon(2),1:2) = AllEle_global(nfcon(2),3:4);
            else
                Tipcoordinate(nfcon(2),1:2) = point1;
            end
            TipStates(nfcon(2)) = nTip+1;
            Index(nfcon(2)) = nAct+1;
            IndexInv(nAct+1) = nfcon(2);
            XfF_global(nfcon(2),:) = Xfinit;
            nAct = nAct + 1;
        else
            if Index(nfcon(2)) < 0.1
                Index(nfcon(2)) = nAct+1;
                IndexInv(nAct+1) = nfcon(2);
                nAct = nAct + 1;
            end
        end
        
        
      %  isCross = 0;
        if isCross > 0.1
            CritThetaa = initAngle;
            disp('******************************************');
            disp(' Fracture Cross natural fracture');
            disp('******************************************');
            TipStates(nf) = TipStates(nf)-999;
            Dnfinit = 0;
            DynamicGrowth_2(1,Pinit,Dnfinit,nf,CritThetaa);
            AllEle_global(nAllEle_global,10) = 1;
            nTip= nTip + 1;
        end
        if isCross < 0.1
            TipStates(nf) = -11;
        end
        HasInter = HasInter-1;
    end
    
end
nActTip = nTip;
clear AllEle ConnList Pres DD;
clear ConnList1 ConnList2 ConnList3;
end