function DynamicGrowth(iscross, Pinit,Dinit,nf,theta)
global ArrivalNF;
global stepL Mat;
global GrowthN Tipcoordinate;
global HasInter MaxEle;
global TipStates Hslurry_global isMechActive_global Hbanking_global CpF_global XfF_global;
global IndexInv nAct Index nAllEle_global ConnList_global AllEle_global DD_global PresF_global;
Pres = zeros(nAct,1);
hslurry = 0*Hslurry_global(1:nAct,:);
hbanking = 0*Hbanking_global(1:nAct,:);
Cp = 0*CpF_global(1:nAct,:);
Xf = 0*XfF_global(1:nAct,:);
DD = zeros(nAct*2,1);
Dn = zeros(nAct,1);
nAct0 = nAct;
for i = 1 : nAct
    Pres(i) = PresF_global(IndexInv(i));
    hslurry(i,:) = Hslurry_global(IndexInv(i),:);
    hbanking(i,:) = Hbanking_global(IndexInv(i),:);
    Cp(i,:) = CpF_global(IndexInv(i),:);
    Xf(i,:) = XfF_global(IndexInv(i),:);
    DD(i) = DD_global(IndexInv(i),:);
    DD(i+nAct) = DD_global(IndexInv(i)+MaxEle);
    Dn(i) = DD_global(IndexInv(i)+MaxEle);
end

epsDis = stepL*1.5;
nAllEle = nAllEle_global;
ConnList = ConnList_global(1:nAllEle_global,:);
AllEle = AllEle_global(1:nAllEle_global,:);
nAllElenew = nAllEle + GrowthN;
AllElenew = [AllEle;zeros(GrowthN,11)];

[~,nc] = size(ConnList);
ConnListnew  = [ConnList;zeros(GrowthN,nc)];
TipStates(nAllEle+GrowthN) = TipStates(nf);
state0 = TipStates(nf);
TipStates(nf) = 0;
Pnew = Pres;
count = 0;
Tipele = nf;
mindis = 1e6;
addL = stepL;
if ConnList(Tipele,2) < 2.2
    if AllEle(nf,10) < 1.1
        [~,~,~,~,mindis1] = CheckCross(Tipele,AllEle(Tipele,1:4),nAllEle,AllEle,ConnList,epsDis);
        [~,~,~,~,mindis2] = CheckCross_HF(Tipele,AllEle(Tipele,1:4),nAllEle,AllEle,ConnList,epsDis,nAct0);
        mindis = min(mindis1,mindis2);
    end
    if mindis < epsDis
        addL = mindis;%mindis/2;
    end
end
connectNF = 0;
left = ConnListnew(Tipele,3);
right = ConnListnew(Tipele,4);
if left > -0.1
    if AllEle(left,10) > 1.1
        connectNF = 1;
    end
end
if right > -0.1
    if AllEle(right,10) > 1.1
        connectNF = 1;
    end
end
ConnListnew(nAllEle+1:nAllEle+GrowthN,1) = ConnList(Tipele,1);
count = count + GrowthN;
% addL = AllEle(Tipele,7)/2;
iright = ConnList(Tipele,4);
ileft =  ConnList(Tipele,3);
TipCoord = Tipcoordinate(Tipele,:);
nfcon = zeros(2,1);
if ConnList(Tipele,2) < 2.1
    if iright < 0.1 || Index(iright) < 0.1
        ConnListnew(Tipele,4) = nAllEle+1;
    else
        ConnListnew(Tipele,3) = nAllEle+1;
    end
else
    countp = 1;
    for i = 3: ConnList(Tipele,2)+2
        BndEle = ConnList(Tipele,i);
        hasfind = 0;
        %         if isNF < 0.1
        if BndEle < 0.1 % && Index(BndEle) < 0.1
            %if ConnList(BndEle,1) == ConnList(Tipele,1)
            %ConnListnew(Tipele,2) = ConnList(Tipele,2)+1;
            hasfind = 1;
            ConnListnew(Tipele,i) = nAllEle+1;
            % Coord1 = AllEle(Tipele,1:2);
            % isin =  isDotIn(Coord1,AllEle(ConnList(Tipele,i),1:4));
            % if isin > 0.1
            ConnListnew(nAllEle+1:nAllEle+GrowthN,1) = ConnList(Tipele,1);
            % else
            %     ConnListnew(nAllEle+1:nAllEle+GrowthN,1) = ConnList(ConnList(Tipele,4),1);
            % end
        else
            %             if AllEle(BndEle,10) < 2.9
            %                 TipStates(BndEle) = -2;
            %             end
            nfcon(countp) = BndEle;
            countp = countp +1 ;
        end
    end
    ConnListnew(Tipele,2)= ConnList(Tipele,2) + 1;
    if hasfind < 0.1
        ConnListnew(Tipele,ConnListnew(Tipele,2)+2) = nAllEle+1;
        ConnListnew(nAllEle+1:nAllEle+GrowthN,1) = ConnList(Tipele,1);
    end
    %end
    %         end
end


alpha = asind(AllEle(Tipele,5));
beta = alpha + theta;
if beta > 360
    beta = beta - 360;
end


%%
Dirvec = [cosd(beta),sind(beta)];

AllElenew(nAllEle+1:nAllEle+GrowthN,7) = addL;

ConnListnew(nAllEle+1:nAllEle+GrowthN,2) = 2;
ConnListnew(nAllEle+1,3) = Tipele;
ConnListnew(nAllEle+1,4) = nAllEle+2;
if ConnList(Tipele,2) > 2.1
    isIn = isDotIn(TipCoord,AllEle(nfcon(1),1:4));
    if isIn > 0.1
        ConnListnew(nAllEle+1,2) = 4;
        ConnListnew(nAllEle+1,5) = nfcon(1);
        ConnListnew(nAllEle+1,6) = nfcon(2);
    end
end

AllElenew(nAllEle+1,1:2) = TipCoord;
AllElenew(nAllEle+1,3:4) = TipCoord+Dirvec*addL;
for i = 2 : GrowthN
    ConnListnew(nAllEle+i,3) = nAllEle+i-1;
    ConnListnew(nAllEle+i,4) = nAllEle+i+1;
    AllElenew(nAllEle+i,1:2) = TipCoord+Dirvec*addL*(i-1);
    AllElenew(nAllEle+i,3:4) = TipCoord+Dirvec*addL*i;
end
if GrowthN < 1.5
    i = 1;
end
ConnListnew(nAllEle+i,4) = -1;
Tipcoordinate(nAllEle+GrowthN,:) = TipCoord + Dirvec*addL*GrowthN;

if beta < -90
    beta = beta + 180;
    for i = 1 : GrowthN
        temp = ConnListnew(nAllEle+i,3);
        ConnListnew(nAllEle+i,3) = ConnListnew(nAllEle+i,4);
        ConnListnew(nAllEle+i,4) = temp;
        temp = AllElenew(nAllEle+i,1:2);
        AllElenew(nAllEle+i,1:2) = AllElenew(nAllEle+i,3:4);
        AllElenew(nAllEle+i,3:4) = temp;
        
    end
else
    if beta > 90  && beta <= 270
        beta = beta - 180;
        for i = 1 : GrowthN
            temp = ConnListnew(nAllEle+i,3);
            ConnListnew(nAllEle+i,3) = ConnListnew(nAllEle+i,4);
            ConnListnew(nAllEle+i,4) = temp;
            temp = AllElenew(nAllEle+i,1:2);
            AllElenew(nAllEle+i,1:2) = AllElenew(nAllEle+i,3:4);
            AllElenew(nAllEle+i,3:4) = temp;
        end
    end
    if beta > 270
        beta = beta - 360;
    end
end
%% Prevent strong curving
%calc= 0;
if iscross < 0.1 && ConnListnew(Tipele,2) < 2.1 %&& calc > 0.1%&& connectNF < 0.1
    epsBeta = 8;
    beta0 = asind(AllElenew(Tipele,5));
    ischange = 0;
    alpha = 0.8;
    if beta0*beta >= 0
        if abs(beta0-beta) > epsBeta
            beta = alpha*beta0+(1-alpha)*beta;
            % Change Original
            ischange = 1;
        end
    else
        if (abs(beta0)+abs(beta) < 180 - epsBeta && abs(beta0)+abs(beta) > 92)
            delta = 180 - (abs(beta0)+abs(beta));
            delta = delta*(1-alpha);
            beta0_2 = beta0 + sign(beta0)*delta;
            if abs(beta0_2) > 90
                beta = -sign(beta0)*(90 - (abs(beta0_2)-90));
            else
                beta = beta0_2;
            end
            ischange = 1;
        end
        if (abs(beta0)+abs(beta) > epsBeta && abs(beta0)+abs(beta) < 20)
            ischange = 1;
            delta = abs(beta0)+abs(beta);
            delta =delta*(1-alpha);
            beta = beta - sign(beta0)*delta/2;
        end
    end
    %ischange = 0;
    if ischange > 0.1
        L = stepL*0.5;
        point1 = AllElenew(Tipele,1:2);
        point2 = AllElenew(Tipele,3:4);
        L0 = AllElenew(Tipele,7);
        dir = [cosd(beta),sind(beta)];
        if abs(point1(1)-TipCoord(1)) < 1e-10 && abs(point1(2)-TipCoord(2)) < 1e-10
            vec1 = -[TipCoord(1)-point2(1),TipCoord(2)-point2(2)];
            if dot(dir,vec1) > 1e-16
                AllElenew(Tipele,1:2) = AllElenew(Tipele,3:4) - L0*dir;
                AllElenew(nAllEle+GrowthN,3:4) = AllElenew(Tipele,3:4) - L0*dir;
                AllElenew(nAllEle+GrowthN,1:2) = AllElenew(Tipele,3:4) -  (L0+L)*dir;
                Tipcoordinate(nAllEle+GrowthN,:) = AllElenew(nAllEle+GrowthN,1:2);
            else
                AllElenew(Tipele,1:2) = AllElenew(Tipele,3:4) + L0*dir;
                AllElenew(nAllEle+GrowthN,1:2) = AllElenew(Tipele,3:4) + L0*dir;
                AllElenew(nAllEle+GrowthN,3:4) = AllElenew(Tipele,3:4) +  (L0+L)*dir;
                Tipcoordinate(nAllEle+GrowthN,:) = AllElenew(nAllEle+GrowthN,3:4);
            end
        else
            vec1 = -[TipCoord(1)-point1(1),TipCoord(2)-point1(2)];
            if dot(dir,vec1) > 1e-16
                AllElenew(Tipele,3:4) = AllElenew(Tipele,1:2) - L0*dir;
                AllElenew(nAllEle+GrowthN,3:4) = AllElenew(Tipele,1:2)  - L0*dir;
                AllElenew(nAllEle+GrowthN,1:2) = AllElenew(Tipele,1:2)  -  (L0+L)*dir;
                Tipcoordinate(nAllEle+GrowthN,:) = AllElenew(nAllEle+GrowthN,1:2);
            else
                AllElenew(Tipele,3:4) =AllElenew(Tipele,1:2)  + L0*dir;
                AllElenew(nAllEle+GrowthN,1:2) = AllElenew(Tipele,1:2)  + L0*dir;
                AllElenew(nAllEle+GrowthN,3:4) = AllElenew(Tipele,1:2)  + (L0+L)*dir;
                Tipcoordinate(nAllEle+GrowthN,:) = AllElenew(nAllEle+GrowthN,3:4);
            end
        end
        AllElenew(nAllEle+GrowthN,7) = L;
        AllElenew(Tipele,5) = sind(beta);
        AllElenew(Tipele,6) = cosd(beta);
        AllElenew(Tipele,8) = AllElenew(Tipele,1)*0.5+AllElenew(Tipele,3)*0.5;
        AllElenew(Tipele,9) = AllElenew(Tipele,2)*0.5+AllElenew(Tipele,4)*0.5;
       % AllElenew(nAllEle+i,9) = AllElenew(Tipele,2)*0.5+AllElenew(Tipele,4)*0.5;
    end
    %%
end

for i = 1 : GrowthN
    AllElenew(nAllEle+i,8) = AllElenew(nAllEle+i,1)*0.5+AllElenew(nAllEle+i,3)*0.5;
    AllElenew(nAllEle+i,9) = AllElenew(nAllEle+i,2)*0.5+AllElenew(nAllEle+i,4)*0.5;
    AllElenew(nAllEle+i,5) = sind(beta);
    AllElenew(nAllEle+i,6) = cosd(beta);
    AllElenew(nAllEle+i,10) = 1;
end
if count ~= GrowthN
    error('Count Error!');
end

%与天然裂缝的相交
%按照Tip Square的分布赋予初始开度
Dshf = DD(Index(Tipele));
Dnhf = DD(Index(Tipele)+nAct);
Dadd = zeros(nAllElenew-nAllEle,2);
xtip = GrowthN-0.5:-1:0.5;
xtip = xtip'*addL;
Li = GrowthN*addL+AllEle(Tipele,7)/2;
Dadd(:,1) = Dshf/sqrt(Li)*sqrt(xtip);
Dadd(:,2) = Dnhf/sqrt(Li)*sqrt(xtip);
Dnew = [DD(1:nAct);Dadd(:,1)*0;DD(1+nAct:nAct*2);Dadd(:,2)*0+Dinit];

Pnew = [Pnew;zeros(GrowthN,1)+Pinit];
hbanking = [hbanking;0*hbanking(1:GrowthN,:)];
hslurry = [hslurry; 0*hslurry(1:GrowthN,:)+Mat.h];
for ii = 1 : GrowthN
    Cp = [Cp;0*Cp(ii,:)+1e-5*Cp(Index(Tipele),:)];
    Xf = [Xf;0*Xf(ii,:)+Xf(Index(Tipele),:)];
end
nAct = nAct + GrowthN;
Index(nAllElenew-GrowthN+1:nAllElenew) = nAct-GrowthN+1:nAct;
IndexInv(nAct-GrowthN+1:nAct) = nAllElenew-GrowthN+1:nAllElenew;
% Check Cross with NF
addnf = 0;
if AllEle(nf,10) < 1.1 && ConnList(Tipele,2) < 2.2
    if nAct > 98
        f = 1;
    end
    [iNF,jNF,Point,isMerge,mindis] = CheckCross(nAllEle+GrowthN,AllElenew(nAllEle+GrowthN,1:4),nAllEle,AllElenew,ConnListnew,epsDis);
    if iNF < -0.1
        iNF = jNF;
    end
    if jNF < -0.1
        jNF = iNF;
    end
    if isMerge > 0.1
        disp('Have Merge');
        HasInter  = HasInter+1;
        if iNF ~= jNF
            CrossType = 2;
            if jNF == ConnList(iNF,3)
                temp = iNF;
                iNF = jNF;
                jNF = temp;
            end
        else
            CrossType = 1;
        end
        % Refine average length
       % addL = stepL/3;
        addL = stepL/2;%stepL/3;
        nAdd = ceil(mindis/addL);
        addL = mindis/nAdd;
        thetai = linspace(-pi,-pi/2,nAdd+1);
        Lseq = zeros(nAdd,1);
        for iL = 1 : nAdd
            Lseq(nAdd-iL+1) = (cos(thetai(iL+1)) + 1)*mindis -(cos(thetai(iL)) + 1)*mindis;
        end
        Dadd = zeros(nAdd,2);
        
        Lseq = ones(nAdd,1)*addL;
        [nAllElenew,AllElenew,ConnListnew] = BuildNewEle(Point,nAllElenew,AllElenew,ConnListnew,CrossType,iNF,jNF,Lseq);
        
        le = ConnListnew(nAllEle+GrowthN,3);
        if le < 0.1
            ConnListnew(nAllEle+GrowthN,3) = nAllEle+GrowthN+1;
        else
            ConnListnew(nAllEle+GrowthN,4) = nAllEle+GrowthN+1;
        end
        iparent = ConnList(iNF,1);
        ArrivalNF(iparent)=1;
        if iNF ~= jNF
            Dnew = [Dnew(1:nAct);Dadd(:,1);
                Dnew(nAct+1:nAct*2);Dadd(:,2)];
            Pnew = [Pnew(1:nAct);zeros(nAdd,1)+Pinit];
        else
            Dnew = [Dnew(1:nAct);Dadd(:,1);
                Dnew(nAct+1:nAct*2);Dadd(:,2)];
            Pnew = [Pnew(1:nAct);zeros(nAdd,1)+Pinit];
        end
        hbanking = [hbanking;0*hbanking(1:nAdd,:)];
        hslurry = [hslurry; 0*hslurry(1:nAdd,:)+Mat.h];
        for ii = 1 : nAdd
            Cp = [Cp;0*Cp(ii,:)+1e-5*Cp(Index(Tipele),:)];
            Xf = [Xf;0*Xf(ii,:)+Xf(Index(Tipele),:)];
        end
        TipStates(nAllEle+GrowthN + nAdd)= 999+state0 ;
        TipStates(nAllEle+GrowthN) = 0;
        nAct = nAct + nAdd;
        Index(nAllEle+GrowthN+1:nAllEle+GrowthN+nAdd) = nAct-nAdd+1:nAct;
        IndexInv(nAct-nAdd+1:nAct) = nAllEle+GrowthN+1:nAllEle+GrowthN + nAdd;
        Tipcoordinate(nAllEle+GrowthN+nAdd,:) = Point;
        if Index(iNF) > 0.1
            addnf = nAllElenew - (nAllEle+GrowthN+nAdd);
            for inf = 1 : addnf
                CpF_global(nAllEle+GrowthN+nAdd+inf,:) = Cp(Index(Tipele),:);
                XfF_global(nAllEle+GrowthN+nAdd+inf,:) = Xf(Index(Tipele),:);
                Hbanking_global(nAllEle+GrowthN+nAdd+inf,:) = Hbanking_global(iNF,:);
                Hslurry_global(nAllEle+GrowthN+nAdd+inf,:) = Hslurry_global(iNF,:);
            end
            nAct = nAct+ addnf;
            Index(nAllEle+GrowthN+nAdd+1:nAllEle+GrowthN+nAdd+addnf) = nAct+1-addnf:nAct;
            IndexInv(nAct+1-addnf:nAct) = nAllEle+GrowthN+nAdd+1:nAllEle+GrowthN+nAdd+addnf;
            isMechActive_global(nAllEle+GrowthN+nAdd+1:nAllEle+GrowthN+nAdd+addnf) = isMechActive_global(iNF);
        end
    end
      % Check Cross with HF
    [iHF,jHF,Point,isMerge,mindis] = CheckCross_HF(nAllEle+GrowthN,AllElenew(nAllEle+GrowthN,1:4),nAllEle,AllElenew,ConnListnew,epsDis,nAct0);
    if isMerge > 0.1
        if iHF ~= jHF
            CrossType = 2;
            if jHF == ConnList(iHF,3)
                temp = iHF;
                iHF = jHF;
                jHF = temp;
            end
        else
            CrossType = 1;
        end
        addL = stepL;%stepL/3;
        nAdd = ceil(mindis/addL);
        addL = mindis/nAdd;
        thetai = linspace(-pi,-pi/2,nAdd+1);
        Lseq = zeros(nAdd,1);
        for iL = 1 : nAdd
            Lseq(nAdd-iL+1) = (cos(thetai(iL+1)) + 1)*mindis -(cos(thetai(iL)) + 1)*mindis;
        end
        Lseq = ones(nAdd,1)*addL;
        
        [nAllElenew,AllElenew,ConnListnew] = BuildNewEle(Point,nAllElenew,AllElenew,ConnListnew,CrossType,iHF,jHF,Lseq);
        
%         le = ConnListnew(nAllEle+GrowthN,3);
%         if le < 0.1
%             ConnListnew(nAllEle+GrowthN,3) = nAllEle+GrowthN+1;
%         else
%             ConnListnew(nAllEle+GrowthN,4) = nAllEle+GrowthN+1;
%         end
%         
        Dadd = zeros(nAdd,2);
        Dadd(:,1) = DD(Index(iHF));
        Dadd(:,2) = Dn(Index(iHF));
        Padd = Pres(Index(iHF));
        
        Dnew = [Dnew(1:nAct);Dadd(:,1);
            Dnew(nAct+1:nAct*2);Dadd(:,2)];
        Pnew = [Pnew(1:nAct);zeros(GrowthN,1)+Padd];
        
        hbanking = [hbanking;0*hbanking(1:GrowthN,:)];
        hslurry = [hslurry; 0*hslurry(1:GrowthN,:)+Mat.h];
        for ii = 1 : GrowthN
            Cp = [Cp;0*Cp(ii,:)+1e-5*Cp(Index(Tipele),:)];
            Xf = [Xf;0*Xf(ii,:)+Xf(Index(Tipele),:)];
        end
        TipStates(nAllEle+GrowthN) = 0;
        nAct = nAct + GrowthN;
        Index(nAllEle+GrowthN+1:nAllEle+GrowthN+nAdd) = nAct-GrowthN+1:nAct;
        IndexInv(nAct-GrowthN+1:nAct) = nAllEle+GrowthN+1:nAllEle+GrowthN+nAdd;
    end
end
AllEle_global(1:nAllElenew,:) = AllElenew;
nAllEle_global = nAllElenew;
ConnList_global(1:nAllElenew,:)  = ConnListnew;
for i = 1 : nAct-addnf
    CpF_global(IndexInv(i),:) = Cp(i,:);
    XfF_global(IndexInv(i),:) = Xf(i,:);
    Hbanking_global(IndexInv(i),:) = hbanking(i,:);
    Hslurry_global(IndexInv(i),:) = hslurry(i,:);
    DD_global(IndexInv(i)) = Dnew(i);
    DD_global(IndexInv(i)+MaxEle) = Dnew(i+nAct-addnf);
    PresF_global(IndexInv(i)) = Pnew(i);
end
end