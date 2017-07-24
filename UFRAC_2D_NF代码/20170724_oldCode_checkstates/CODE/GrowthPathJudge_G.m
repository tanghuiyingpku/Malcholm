function   GrowthPathJudge_G(dact,Dnfinit)
disp('** Judge the intersection behavior with Energy Release Rate Criteria')
global KIChf KICnf
global TipStates;
global nTip;
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
        end
        [Gnf,num] = max(G);
        Gnf = Gnf/KICnf;
        % HF
        Pinit = Mat.Pp;
        Xfinit = Xf(i,:);
        Ghf = (K1^2 + K2^2)/Mat.E/KIChf;
        if Gnf > Ghf  && Gnf > 1.01
            HasInter = 0;
            disp('******************************************');
            disp(' Fracture turn into natural fracture');
            disp('******************************************');
            if abs(G(1) - G(2)) < 1E-6
                disp(' Enter Both Sides');
                AllEle(nfcon(num),10) = 3;
                AllEle(nfcon(3-num),10)  = 3;
                TipStates(nfcon(num)) = nTip+1;
                TipStates(nfcon(3-num)) = nTip+2;
                nTip = nTip + 2;
                Pnew(nAct+1:nAct+2) = Pinit;
                Dn(nAct+1:nAct+2) = Dnfinit;
                Xf(nAct+1:nAct+2,:) = Xfinit;
                Index(nfcon(num)) = nAct+1;
                Index(nfcon(3-num)) = nAct+2;
                IndexInv(nAct+1) = nfcon(num);
                IndexInv(nAct+2) = nfcon(3-num);
                nAct = nAct + 2;
            else
                disp(' Only Enter one Sides');
                AllEle(nfcon(num),10) = 3;
                left = ConnList(nfcon(num),3);
                right = ConnList(nfcon(num),4);
                if left ~= nfcon(3-num)
                    AllEle(left,10) = 3;
                    TipStates(left) =  nTip + 1;
                    TipStates(right) = -2;
                    Pnew(nAct+1:nAct+2) = Pinit;
                    Dn(nAct+1:nAct+2) = Dnfinit;
                    Xf(nAct+1,:) = Xfinit;
                    Xf(nAct+2,:) = Xfinit;
                    Index(nfcon(num)) = nAct+1;
                    Index(left) = nAct+2;
                    IndexInv(nAct+1) = nfcon(num);
                    IndexInv(nAct+2) = left;
                    nAct = nAct + 2;
                else
                    AllEle(right,10) = 3;
                    TipStates(right) =  nTip + 1;
                    TipStates(left) = -2;
                    Pnew(nAct+1:nAct+2) = Pinit;
                    Dn(nAct+1:nAct+2) = Dnfinit;
                    Xf(nAct+1,:) = Xfinit;
                    Xf(nAct+2,:) = Xfinit;
                    Index(nfcon(num)) = nAct+1;
                    Index(right) = nAct+2;
                    IndexInv(nAct+1) = nfcon(num);
                    IndexInv(nAct+2) = right;
                    nAct = nAct + 2;
                end
                nTip = nTip + 1;
            end
            TipStates(IndexInv(i)) = -1;
            DD = [Ds(1:nAct);Dn(1:nAct)];
            Pres = Pnew(1:nAct);
        else
            if Ghf > 1.01
                HasInter = -1;
                [CritTheta,~] = FindCritDirSingle(K1,K2,type);
                disp('******************************************');
                disp(' Fracture Cross natural fracture');
                disp('******************************************');
                TipStates(i) = TipStates(i)-999;
                Dinit = Dnfinit;
                DynamicGrowth(Pinit,Dinit,IndexInv(i),CritTheta);
                % DynamicGrowth(Mat,i,Fractures,nAllEle,ConnList,AllEle,CritTheta,Mat.Sxx,Mat.Syy,Mat.Sxy,DD(1:nAllEle*2,:),Pres(1:nAllEle,1));
                TipStates(IndexInv(i)) = -1;
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
