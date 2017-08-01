function [RelaK,dt0] = GrowthPathJudge_VirtualElemnet(CritTheta,dt0,CurT)
disp('** Judge the intersection behavior with Virtual Element Method')
% global variables
global  Tipcoordinate TipStatesInv  Mat KICnf XfF_global
global TipStates ;
global nTip 
global IndexInv nAct Index MaxEle
global DD_global PresF_global AllEle_global ConnList_global
global HasInter
% Temperary Values
nEle = nAct + 1; % Active one element each time
DD = zeros(nEle*2,1);
Pres = zeros(nEle,1);
AllEle = zeros(nEle,11);
ConnList = zeros(nEle,8);
leftT = -1e6;
rightT = 1e6;
% Set values for the existing elements
for i = 1 : nAct
    DD(i) = DD_global(IndexInv(i));
    DD(i+nEle) = DD_global(MaxEle+IndexInv(i));
    Pres(i) = PresF_global(IndexInv(i));
    AllEle(i,:)  = AllEle_global(IndexInv(i),:);
    ConnList(i,:) = ConnList_global(IndexInv(i),:);
end
% Three independent path
ConnList1 = ConnList;
ConnList2 = ConnList;
ConnList3 = ConnList;

% [DD0,P0] = GetInitialDP();%calcEnergy();
% Gn0 = Gn0 - G_inj;
% Find possible connecting elements and try the virtual element - calc
% opening based on given pressure then compare the relative stress
% intensity factor
for i = 1 : nTip
    nf = TipStatesInv(i);
    if TipStates(nf) > 998
        [Tiptype,~] = TipType(nf,AllEle_global,ConnList_global);
        type0 = ConnList_global(nf,1);
        con = ConnList_global(nf,3:5);
        nfcon = zeros(2,1);
        count = 1;
        cpoint = [AllEle_global(nf,3),AllEle_global(nf,4)];
        if Tiptype > 1.1
            cpoint = [AllEle_global(nf,1),AllEle_global(nf,2)];
        end
        for j = 1 : 3
            typej = ConnList_global(con(j),1);
            if typej ~= type0
                nfcon(count) = con(j);%global numbering
                count = count + 1;
            end
        end
        Pinit = PresF_global(nf)+ 0;
        Xfinit = XfF_global(nf,:);
        initAngle = CritTheta(i);
        MeetCriteria = 0;
        while MeetCriteria < 0.1
            RelaK = zeros(3,1);
            %% NF1
            Pres(nEle) = Pinit;
            AllEle(nEle,:) = AllEle_global(nfcon(1),:);
            d =  AllEle(Index(nf),7);
            Tipcoord = Tipcoordinate(nf,:);
            ta = asind(AllEle(nEle,5));
            if abs(Tipcoord(1) - AllEle(nEle,1)) < 1e-5
                AllEle(nEle,3) = AllEle(nEle,1) + d*cosd(ta);
                AllEle(nEle,4) = AllEle(nEle,2) + d*sind(ta);
            else
                AllEle(nEle,1) = AllEle(nEle,3) - d*cosd(ta);
                AllEle(nEle,2) = AllEle(nEle,4) - d*sind(ta);
            end
            AllEle(nEle,8) = 0.5*(AllEle(nEle,1) + AllEle(nEle,3));
            AllEle(nEle,9) = 0.5*(AllEle(nEle,2) + AllEle(nEle,4));
            ConnList1(nEle,:) = ConnList_global(nfcon(1),:);
            nconn = ConnList_global(nf,2);
            ihf = 0;
            for jj =  1: nconn
                nb = ConnList_global(nf,2+jj);
                if AllEle_global(nb,10) < 1.1
                    ihf = nb;
                end
            end
            ConnList1(Index(nf),2) =2;
            ConnList1(Index(nf),3) = ihf;
            ConnList1(Index(nf),4) = nfcon(1);
               
            disp('********************Test NF branch 1********************')
            DD1 = CalcDD_Cross(dt0,CurT,nEle,DD,Pres,AllEle,ConnList1);
            %DrawDisplacement_local(nEle,DD1,AllEle,ConnList1,Tipcoord,d/2);
            % Get KIC at this nf element
            Dn = DD1(nEle*2);
            Ds = DD1(nEle);
            disp('********************Test NF branch 2********************')
            KI11 = -Mat.G*1e3/4/(1-Mat.miu)*(2*pi/d)^0.5*Dn;
            KI21 =abs( -Mat.G*1e3/4/(1-Mat.miu)*(2*pi/d)^0.5*Ds);
            RelaK(1) = KI11/Mat.K1NF;%(KI11/Mat.K1NF)^2 + (KI21/Mat.K2NF)^2;
            
            %% Virtual Element NF(2)
            Pres(nEle) = Pinit;
            AllEle(nEle,:) = AllEle_global(nfcon(2),:);
            % d =  AllEle_global(nfcon(2),7);
            AllEle(nEle,7) = d;
            Tipcoord = Tipcoordinate(nf,:);
            ta = asind(AllEle(nEle,5));
            if abs(Tipcoord(1) - AllEle(nEle,1)) < 1e-5
                AllEle(nEle,3) = AllEle(nEle,1) + d*cosd(ta);
                AllEle(nEle,4) = AllEle(nEle,2) + d*sind(ta);
            else
                AllEle(nEle,1) = AllEle(nEle,3) - d*cosd(ta);
                AllEle(nEle,2) = AllEle(nEle,4) - d*sind(ta);
            end
            AllEle(nEle,8) = 0.5*(AllEle(nEle,1) + AllEle(nEle,3));
            AllEle(nEle,9) = 0.5*(AllEle(nEle,2) + AllEle(nEle,4));
            ConnList2(nEle,:) = ConnList_global(nfcon(2),:);
            ConnList1(Index(nf),2) =2;
            ConnList1(Index(nf),3) = ihf;
            ConnList1(Index(nf),4) = nfcon(2);
            DD2 = CalcDD_Cross(dt0,CurT,nEle,DD,Pres,AllEle,ConnList2);
           % DrawDisplacement_local(nEle,DD2,AllEle,ConnList2,Tipcoord,d/2);
            % Get KIC at this nf element
            Dn = DD2(nEle*2);
            Ds = DD2(nEle);
            disp('********************Test NF branch 2********************')
            KI12 = -Mat.G*1e3/4/(1-Mat.miu)*(2*pi/d)^0.5*Dn;
            KI22 =abs( -Mat.G*1e3/4/(1-Mat.miu)*(2*pi/d)^0.5*Ds);
            RelaK(2) = KI12/Mat.K1NF;%(KI12/Mat.K1NF)^2 + (KI22/Mat.K2NF)^2;
           %% HF
            Pres(nEle) = Pinit;
            AllEle(nEle,:) = AllEle(Index(nf),:);
            ConnList(nEle,:) = ConnList(Index(nf),:);
            %d = AllEle_global(nf,7);
            alpha = asind(AllEle(Index(nf),5));
            beta = alpha + initAngle;
            if beta > 360
                beta = beta - 360;
            end
            Dirvec = [cosd(beta),sind(beta)];
            TipCoord = Tipcoordinate(nf,:);
            AllEle(nEle,1:2) = TipCoord;
            AllEle(nEle,3:4) = TipCoord+Dirvec*d;
            if beta < -90
                beta = beta + 180;
                temp = AllEle(nEle,1:2);
                AllEle(nEle,1:2) = AllEle(nEle,3:4);
                AllEle(nEle,3:4) = temp;
            else
                if beta > 90  && beta <= 270
                    beta = beta - 180;
                    temp = AllEle(nEle,1:2);
                    AllEle(nEle,1:2) = AllEle(nEle,3:4);
                    AllEle(nEle,3:4) = temp;
                end
                if beta > 270
                    beta = beta - 360;
                end
            end
            AllEle(nEle,5) = sind(beta);
            AllEle(nEle,6) = cosd(beta);
            AllEle(nEle,8) = 0.5*(AllEle(nEle,1) + AllEle(nEle,3));
            AllEle(nEle,9) = 0.5*(AllEle(nEle,2) + AllEle(nEle,4));
            ConnList3(Index(nf),2) =2;
            if  ConnList3(Index(nf),3) < 0.1
                ConnList3(Index(nf),3) =MaxEle;
            else
                ConnList3(Index(nf),4) =MaxEle;
            end
            Index(MaxEle) = nEle;
            ConnList3(nEle,:) = ConnList(Index(nf),:);
            ConnList3(nEle,2) = 1;
            ConnList3(nEle,3) = nf;
            ConnList3(nEle,4) = -1;
            DD3 = CalcDD_Cross(dt0,CurT,nEle,DD,Pres,AllEle,ConnList3);
            % Get KIC at this nf element
            Dn = DD3(nEle*2);
            Ds = DD3(nEle);
            KI1 = -Mat.G*1e3/4/(1-Mat.miu)*(2*pi/d)^0.5*Dn;
            KI2 = -Mat.G*1e3/4/(1-Mat.miu)*(2*pi/d)^0.5*Ds;
            a_theta = -90;
            b_theta = 90;
            [K3,~]= FindMaxTheta(a_theta,b_theta,KI1,KI2);
            RelaK(3) = K3;
%             %
%             if RelaK(1) > RelaK(2)
%                 a_theta = -90;
%                 b_theta = 90;
%                 [Kmax,~]= FindMaxTheta(a_theta,b_theta,KI11,KI21);
%             else
%                 a_theta = -90;
%                 b_theta = 90;
%                 [Kmax,~]= FindMaxTheta(a_theta,b_theta,KI12,KI22);
%             end
%             RelaK(3) = Kmax;
            Kmax = max(RelaK);
            MeetCriteria = 0;
            if abs(Kmax - 1) <= 0.1
                MeetCriteria = 1;
            else
                [leftT,rightT,dt0]  = Adjustdt(leftT,rightT,dt0,Kmax);
                disp('Change Change Time Step')
            end
        end
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