function RelaK = GrowthPathJudge_GIGII_backup(CritG,CritTheta,K1c,K2c)
disp('** Judge the intersection behavior with Virtual Element Method')
% global variables
global  Tipcoordinate TipStatesInv  Mat  XfF_global
global TipStates ;
global nTip
global IndexInv nAct Index MaxEle
global DD_global PresF_global AllEle_global  ConnList_global
global HasInter

for i = 1 : nTip
    nf = TipStatesInv(i);
    if TipStates(nf) > 998
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
        % 1
        Ele= AllEle_global(nfcon(1),:);
        Tipcoord = Tipcoordinate(nf,:);
        ta = asind(Ele(5));
        if abs(Tipcoord(1) - Ele(1)) < 1e-5
            angle1 = ta;
            theta1 = angle1 - angleH;
        else
            angle1 = -ta;
            theta1 = angle1 - angleH;
        end
        k11 = 0.5*cosd(theta1/2)*(K1*(1+cosd(theta1))-3*K2*sind(theta1));
        k21 = 0.5*cosd(theta1/2)*(K1*sind(theta1)+K2*(3*cosd(theta1)-1));
        RelaK(1) = (k11/Mat.K1NF)^2 + (k21/Mat.K2NF)^2;
        
         %2
        Ele= AllEle_global(nfcon(2),:);
        Tipcoord = Tipcoordinate(nf,:);
        ta = asind(Ele(5));
        if abs(Tipcoord(1) - Ele(1)) < 1e-5
            angle2 = (ta);
            theta2 = angle2 - angleH;
        else
            angle2 = -(ta);
            theta2 = angle2 - angleH;
        end
        k12 = 0.5*cosd(theta2/2)*(K1*(1+cosd(theta2))-3*K2*sind(theta2));
        k22 = 0.5*cosd(theta2/2)*(K1*sind(theta2)+K2*(3*cosd(theta2)-1));
        RelaK(2) = (k12/Mat.K1NF)^2 + (k22/Mat.K2NF)^2;
        
        RelaK(3) = CritG(i);
        
        Kmax = max(RelaK);
        
        %  DrawDisplacement_local(nEle,DD2,AllEle,ConnList2,Tipcoord,d,nfcon(cf))
        
        if abs(RelaK(1) - Kmax) < 1e-6
            disp('Dilate into NF 1');
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
            PresF_global(nfcon(1)) = Pinit;
            DD_global(MaxEle+nfcon(1)) = 0;
            XfF_global(nfcon(1),:) = Xfinit;
            TipStates(nf) = -1;
            TipStates(nfcon(2)) = -3;
        end
        if abs(RelaK(2) - Kmax) < 1e-6
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
            nAct = nAct + 1;
            PresF_global(nfcon(2)) = Pinit;
            DD_global(MaxEle+nfcon(2)) = 0;
            XfF_global(nfcon(2),:) = Xfinit;
            TipStates(nf) = -1;
            TipStates(nfcon(1)) = -3;
        end
        if abs(RelaK(3) - Kmax) < 1e-6
            disp('Cross NF');
            CritThetaa = initAngle;
            disp('******************************************');
            disp(' Fracture Cross natural fracture');
            disp('******************************************');
            TipStates(nf) = TipStates(nf)-999;
            Dnfinit = 0;
            DynamicGrowth(1,Pinit,Dnfinit,nf,CritThetaa);
            TipStates(nfcon(1)) = -3;
            TipStates(nfcon(2))  = -3;
        end
    end
end
HasInter = -1;
clear AllEle ConnList Pres DD;
clear ConnList1 ConnList2 ConnList3;
end