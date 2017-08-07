
function  GrowthPathJudge_GIGII2(nf0,CritTheta,K1c,K2c)
disp('** Judge the intersection behavior with Virtual Element Method')
% global variables
global  Tipcoordinate TipStatesInv  Mat  XfF_global Fractures
global TipStates ;
global nTip isMechActive_global;
global IndexInv nAct Index
global  PresF_global AllEle_global  ConnList_global nAllEle_global
global HasInter  GrowNumber;
AlreadyAct = zeros(nAllEle_global,1);
GrowNumber2 = GrowNumber;
for i = 1 : nTip
    nf = TipStatesInv(i);
    if nf > 0.1 && TipStates(nf) > 998 && nf == nf0
        K1 = K1c(i);
        K2 = K2c(i);
        [Tiptype,~] = TipType(nf,AllEle_global,ConnList_global);
        type0 = ConnList_global(nf,1);
        con = ConnList_global(nf,3:5);
        for ic = 1 : 3
            if con(ic) < -3
                ConnList_global(nf,2+ic) = -con(ic);
                con(ic) =-con(ic);
            end
        end
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
            Fang = Mat.fric;
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
        
        disp('Dilate into NF 1');
        % Make NF type to be activated NF -- type 3
        Index_Frac = ConnList_global(nfcon(1),1);
        % Active
        %         startI = Fractures{Index_Frac}.index(1);
        %从交点往外扩GrowNumber个
        Actnum = zeros(50,1);
        ActTip = zeros(50,1);
        countN = 0;
        countT = 0;
        count = 0;
        jump = 0;
        past = nfcon(2);
        ic = nfcon(1);
        for ie = 1 : ConnList_global(ic,2)
            bc = ConnList_global(ic,2+ie);
            if bc < -3
                ConnList_global(ic,2+ie) = -bc;
            end
        end
        for ie = 1 : ConnList_global(past,2)
            bc = ConnList_global(past,2+ie);
            if bc < -3
                ConnList_global(past,2+ie) = -bc;
            end
        end
        if  Index(ic) < 0.1
            countN = countN + 1;
            Actnum(countN) = ic;
            while  count <  GrowNumber2
                if min(ConnList_global(ic,3:6)) < -0.1 && min(ConnList_global(ic,3:6)) > -2
                    countT = countT+ 1;
                    ActTip(countT) = ic;
                    jump = 1;
                    break;
                end
                for ie = 1 : ConnList_global(ic,2)
                    bc = ConnList_global(ic,2+ie);
                    if bc < -3
                        ConnList_global(ic,2+ie) = -bc;
                        bc = -bc;
                    end
                    if bc > 0.1 && bc ~= past
                        if ConnList_global(bc,1) == Index_Frac
                            count = count + 1;
                            %Activate current
                            if Index(bc) < 0.1  && AlreadyAct(bc) < 0.1
                                countN = countN + 1;
                                Actnum(countN) = bc;
                                AlreadyAct(bc) = 1;
                            end
                            past = ic;
                            ic = bc;
                        else
                            %cross with other natural fractures
                            if Index(bc) < 0.1  && AlreadyAct(bc) < 0.1
                                countN = countN + 1;
                                Actnum(countN) = bc;
                                countT = countT+ 1;
                                ActTip(countT) = bc;
                                AlreadyAct(bc) = 1;
                            end
                        end
                    end
                end
            end
        end
        if jump < 0.1 && countN > 0.1
            countT = countT+ 1;
            ActTip(countT) = ic;
        end
        count = 0;
        jump = 0;
        past = nfcon(1);
        ic = nfcon(2);
        if  Index(ic) < 0.1
            countN = countN + 1;
            Actnum(countN) = ic;
            while  count <  GrowNumber2
                if min(ConnList_global(ic,3:6)) < -0.1 && min(ConnList_global(ic,3:6)) > -2
                    countT = countT+ 1;
                    ActTip(countT) = ic;
                    jump = 1;
                    break;
                end
                for ie = 1 : ConnList_global(ic,2)
                    bc = ConnList_global(ic,2+ie);
                    if bc < -3
                        ConnList_global(ic,2+ie) = -bc;
                        bc = -bc;
                    end
                    if bc > 0.1 && bc ~= past
                        if ConnList_global(bc,1) == Index_Frac
                            count = count + 1;
                            %Activate current
                            if Index(bc) < 0.1  && AlreadyAct(bc) < 0.1
                                countN = countN + 1;
                                Actnum(countN) = bc;
                                AlreadyAct(bc) = 1;
                            end
                            past = ic;
                            ic = bc;
                        else
                            %cross with other natural fractures
                            if Index(bc) < 0.1 && AlreadyAct(bc) < 0.1
                                countN = countN + 1;
                                Actnum(countN) = bc;
                                countT = countT+ 1;
                                ActTip(countT) = bc;
                                AlreadyAct(bc) = 1;
                            end
                        end
                    end
                end
            end
        end
        if jump < 0.1 && countN > 0.1
            countT = countT+ 1;
            ActTip(countT) = ic;
        end
        for in = 1:countN
            nb = Actnum(in);
            AllEle_global(nb,10) = 3;
            PresF_global(nb) = Mat.Pp;
            XfF_global(nb,:) = Xfinit;
            Index(nb) = nAct+1;
            IndexInv(nAct+1) = nb;
            nAct = nAct + 1;
            isMechActive_global(nAct) = -2;
        end
        for in = 1 : countT
            bc = ActTip(in);
            TipStates(bc) = nTip+1;
            Tipcoordinate(bc,1:2) = FindTipCoord(bc);
            nTip = nTip + 1;
        end
        isCross = 0;
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
            nTip = nTip - 1;
        end
        HasInter = HasInter-1;
    end
    
end

clear AllEle ConnList Pres DD AlreadyAct;
clear ConnList1 ConnList2 ConnList3;
end