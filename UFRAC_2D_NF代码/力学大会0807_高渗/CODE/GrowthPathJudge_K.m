function   GrowthPathJudge_K(dact,Dnfinit)
disp('** Judge the intersection behavior according to Tip stress distribution')
global KIChf Fractures;
global TipStates;
global nTip  Tipcoordinate;
global IndexInv nAct Index MaxEle  XfF_global Fluid
global DD_global PresF_global AllEle_global ConnList_global nAllEle_global Mat
global HasInter

Pnew = zeros(nAct+10,1);
Ds = zeros(nAct+10,1);
Dn = zeros(nAct+10,1);
AllEle = AllEle_global(1:nAllEle_global,:);
ConnList = ConnList_global(1:nAllEle_global,:);
Xf = zeros(nAct+10,Fluid.nfluid);
for i = 1 : nAct
    Xf(i,:) = XfF_global(IndexInv(i),:);
    Pnew(i) = PresF_global(IndexInv(i));
    Ds(i) = DD_global(IndexInv(i));
    Dn(i) = DD_global(IndexInv(i)+MaxEle);
end
Pres = Pnew;
DD = [Ds;Dn];
nAllEle = nAllEle_global;
% TipStates = 999 means the intersection location on HF
for i = nAct-dact+1: nAct
    if TipStates(IndexInv(i)) > 998
        [type,~] = TipType(IndexInv(i),AllEle,ConnList);
        Dss= DD(i);
        Dnn= Dn(i);
        d = AllEle(IndexInv(i),7)/2;
        type0 = ConnList(IndexInv(i),1);
        K1 = -Mat.G*1e3/4/(1-Mat.miu)*(2*pi/d)^0.5*Dnn;
        K2 = -Mat.G*1e3/4/(1-Mat.miu)*(2*pi/d)^0.5*Dss;
        a = -90;
        b = -30;
        [hasFind,c]= FindOpenTheta(a,b,K1,K2);
        while hasFind < 0.1
            a = a+60;
            b = b+60;
            [hasFind,c]= FindOpenTheta(a,b,K1,K2);
        end
        if hasFind < 0.1
            %没有找到临界角度
            error('Fail to Find the critical Angle');
        end
        KC = K1* 0.5*cosd(c/2)*(1+cosd(c))-K2*3/2*cosd(c/2)*sind(c);
        if KC < Mat.Kmin/5
            continue;
        else
            con = ConnList(IndexInv(i),3:5);
            nfcon = zeros(2,1);
            count = 1;
            for j = 1 : 3
                typej = ConnList(con(j),1);
                if typej ~= type0
                    nfcon(count) = con(j);
                    count = count + 1;
                end
            end
            sini = AllEle(IndexInv(i),5);
            thetai = asind(sini);
            cpoint = [AllEle(IndexInv(i),3),AllEle(IndexInv(i),4)];
            if type > 1.1
                cpoint = [AllEle(IndexInv(i),1),AllEle(IndexInv(i),2)];
            end
            %cpoint is the tip coordinates
            vecmain = [cosd(thetai),sind(thetai)];
            % vecmain is the vector point outside the fracture
            initAngle = zeros(2,1);
            rc = zeros(2,1);
            ds = zeros(2,1);
            G = zeros(2,1);
            for j = 1 : 2
                nf = nfcon(j);
                sinnf = AllEle(nf,5);
                alpha = asind(sinnf);
                vecnf = [AllEle(nf,8)-cpoint(1),AllEle(nf,9)-cpoint(2)];
                %[Cos Sin]
                veccoord = [AllEle(nf,6) AllEle(nf,5)];
                if type < 1.1
                    theta0 = alpha - thetai;
                else
                    theta0 = alpha - thetai + 180;
                end
                % Whether along the local coordiate
                vecdir = dot(vecnf,veccoord);
                if vecdir < -1e-10
                    theta0 = theta0 + 180;
                end
                k11 = 0.5*cosd(theta0/2)*(K1*(1+cosd(theta0))-3*K2*sind(theta0));
                k12 = 0.5*cosd(theta0/2)*(K1*sind(theta0)+K2*(3*cosd(theta0)-1));
                G(j) = k11^2/Mat.E+k12^2/Mat.E;
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
                [initAngle(j),rc(j)] = Criteria(thetai,Mat.Sxx,Mat.Syy,Mat.Sxy,theta0,K1,K2,To);
                ds(j) = M_C_stability(thetai,Mat.Sxx,Mat.Syy,Mat.Sxy,theta0,K1,K2,rc(j),Fang,So);
            end
            Pinit = abs(Mat.Syy);
            Xfinit = Xf(i,:);
            if type > 1.1
                initAngle = initAngle  - 180;
            end
            if min(ds) < -1e-16
                % Sliding
                iscross = 0 ;
            else
                iscross = 1;
            end
            [~,num] = min(ds);
            if iscross < 0.1
                HasInter = 0;
                disp('******************************************');
                disp(' Fracture turn into natural fracture');
                disp('******************************************');
                if abs(ds(1) - ds(2)) < 1E-12
                    disp(' Enter Both Sides');
                    
                    addi = 1;
                    if Index(nfcon(num)) < 0.1 && G(num) > G(3-num)
                        % Opened
                        AllEle(nfcon(num),10) = 3;
                        point1 = AllEle(nfcon(num),1:2);
                        dis =  sqrt((point1(1) - cpoint(1))^2 + (point1(2) - cpoint(2))^2 );
                        if dis < 1e-6
                            Tipcoordinate(nfcon(num),1:2) = AllEle(nfcon(num),3:4);
                        else
                            Tipcoordinate(nfcon(num),1:2) = point1;
                        end
                        TipStates(nfcon(num)) = nTip+1;
                        nTip = nTip + 1;
                        Pnew(nAct+addi) = Pinit;
                        Dn(nAct+addi) = Dnfinit;
                        Xf(nAct+addi,:) = Xfinit;
                        Index(nfcon(num)) = nAct+addi;
                        IndexInv(nAct+addi) = nfcon(num);
                        nAct = nAct + 1;
                    end
                    if Index(nfcon(3-num)) < 0.1 && G(num) < G(3-num)
                         AllEle(nfcon(3-num),10)  = 3;
                        point1 = AllEle(nfcon(3-num),1:2);
                        dis =  sqrt((point1(1) - cpoint(1))^2 + (point1(2) - cpoint(2))^2 );
                        if dis < 1e-6
                            Tipcoordinate(nfcon(3-num),1:2) = AllEle(nfcon(3-num),3:4);
                        else
                            Tipcoordinate(nfcon(3-num),1:2) = point1;
                        end
                        TipStates(nfcon(3-num)) = nTip+1;
                        nTip = nTip + 1;
                        Pnew(nAct+addi) = Pinit;
                        Dn(nAct+addi) = Dnfinit;
                        Xf(nAct+addi,:) = Xfinit;
                        Index(nfcon(3-num)) = nAct+addi;
                        IndexInv(nAct+addi) = nfcon(3-num);
                        nAct = nAct + 1;
                    end
                else
                    disp(' Only Enter one Sides');
                    AllEle(nfcon(num),10) = 3;
                    TipStates(nfcon(num)) =  nTip + 1;
                    TipStates(nfcon(3-num)) = -2;
                    Pnew(nAct+1) = Pinit;
                    Dn(nAct+1) = Dnfinit;
                    Xf(nAct+1,:) = Xfinit;
                    point1 = AllEle(nfcon(num),1:2);
                    dis =  sqrt((point1(1) - cpoint(1))^2 + (point1(2) - cpoint(2))^2 );
                    if dis < 1e-6
                        Tipcoordinate(nfcon(num),1:2) = AllEle(nfcon(num),3:4);
                    else
                        Tipcoordinate(nfcon(num),1:2) = point1;
                    end
                    Index(nfcon(num)) = nAct+1;
                    IndexInv(nAct+1) = nfcon(num);
                    nAct = nAct + 1;
                    nTip = nTip + 1;
                end
                TipStates(IndexInv(i)) = -1;
                DD = [Ds(1:nAct);Dn(1:nAct)];
                Pres = Pnew(1:nAct);
            else
                HasInter = -1;
                CritTheta = initAngle(num);
                disp('******************************************');
                disp(' Fracture Cross natural fracture');
                disp('******************************************');
                TipStates(IndexInv(i)) = TipStates(IndexInv(i))-999;
                DynamicGrowth(1,Pinit,Dnfinit,IndexInv(i),CritTheta);
                TipStates(nfcon(num)) = -3;
                TipStates(nfcon(3-num))  = -3;
                % DynamicGrowth(Mat,i,Fractures,nAllEle,ConnList,AllEle,CritTheta,Mat.Sxx,Mat.Syy,Mat.Sxy,DD(1:nAllEle*2,:),Pres(1:nAllEle,1));
               % TipStates(IndexInv(i)) = -1;
            end
        end
    end
end
if HasInter < 0.1 && HasInter > -0.01
    for i = 1 : nAct
        PresF_global(IndexInv(i)) = Pres(i);
        DD_global(IndexInv(i)) = DD(i);
        DD_global(IndexInv(i) + MaxEle) = DD(nAct+i);
        XfF_global(IndexInv(i),:) = Xf(nAct,:);
    end
    AllEle_global(1:nAllEle,:) = AllEle;
    nAllEle_global = nAllEle;
    ConnList_global(1:nAllEle,:) = ConnList;
end
    function [hasFind,c]= FindOpenTheta(a,b,KI1,KI2)
        [c,fc] = BiSearch(a,b,KI1,KI2);
        g = Fnc2ndOrder(c,KI1,KI2);
        if g > 0 && abs(fc) < 1e-8
            hasFind = 1;
        else
            hasFind = 0;
            c  = -1;
        end
    end
    function f = DeltFnc(theta,KI1v,KI2v)
        f = cosd(theta/2)*(KI1v*sind(theta) + KI2v*(3*cosd(theta)-1))/1e6;
    end
    function g = Fnc2ndOrder(theta,KI1v,KI2v)
        g = 2*KI1v*cosd(theta/2)*(3*cosd(theta)-1) - KI2v*sind(theta/2)*(9*cosd(theta)+5);
        g = g/1e6;
    end
    function [c,fc] = BiSearch(a,b,KI1,KI2)
        dlt = 1e-6;
        dlt2 = 1e-10;
        while abs(a-b) > dlt
            c = (a+b)/2;
            fb = DeltFnc(b,KI1,KI2);
            fc = DeltFnc(c,KI1,KI2);
            if abs(fc) <dlt2
                break;
            else
                if fc*fb < 0
                    a = c;
                else
                    b = c;
                end
            end
        end
    end
end
