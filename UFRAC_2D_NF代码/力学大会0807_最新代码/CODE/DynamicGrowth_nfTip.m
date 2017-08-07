function DynamicGrowth_nfTip(isNF, Pinit,Dinit,nf,theta)
global ArrivalNF;
global stepL Mat Tipcoordinate;
global GrowthN nHF_global nFracture; 
global HasInter MaxEle Fractures;
global TipStates Hslurry_global Hbanking_global CpF_global XfF_global;
global IndexInv nAct Index nAllEle_global ConnList_global AllEle_global DD_global PresF_global;
Pres = zeros(nAct,1);
hslurry = 0*Hslurry_global(1:nAct,:);
hbanking = 0*Hbanking_global(1:nAct,:);
Cp = 0*CpF_global(1:nAct,:);
Xf = 0*XfF_global(1:nAct,:);
DD = zeros(nAct*2,1);

for i = 1 : nAct
    Pres(i) = PresF_global(IndexInv(i));
    hslurry(i,:) = Hslurry_global(IndexInv(i),:);
    hbanking(i,:) = Hbanking_global(IndexInv(i),:);
    Cp(i,:) = CpF_global(IndexInv(i),:);
    Xf(i,:) = XfF_global(IndexInv(i),:);
    DD(i) = DD_global(IndexInv(i),:);
    DD(i+nAct) = DD_global(IndexInv(i)+MaxEle);
end

epsDis = stepL*GrowthN*2.1;
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
        [~,~,~,~,mindis] = CheckCross(Tipele,AllEle(Tipele,1:4),nAllEle,AllEle,ConnList,epsDis);
    end
    if mindis < epsDis
        addL = mindis;
    end
end

count = count + GrowthN;
% addL = AllEle(Tipele,7)/2;
iright = ConnList(Tipele,4);
ileft =  ConnList(Tipele,3);
nfcon = zeros(2,1);

if ConnList(Tipele,2) < 2.1
    if iright < 0.1 || Index(iright) < 0.1
        ConnListnew(Tipele,4) = nAllEle+1;
        TipCoord = AllEle(Tipele,3:4);
        ConnListnew(nAllEle+1:nAllEle+GrowthN,1) =nFracture + 1;
        nHF_global = nHF_global +1;
        nFracture = nFracture + 1;
        Fractures{nFracture}.Type = 1;
        if isNF > 0.1
            if AllEle(iright,10) < 2.9
                TipStates(iright) = -2;
            end
            nfcon(1) = iright;
        end
    else
        ConnListnew(Tipele,3) = nAllEle+1;
        TipCoord = AllEle(Tipele,1:2);
        ConnListnew(nAllEle+1:nAllEle+GrowthN,1) =nFracture + 1;
        nHF_global = nHF_global +1;
        nFracture = nFracture + 1;
        Fractures{nFracture}.Type = 1;
        if isNF > 0.1
            if AllEle(ileft,10) < 2.9
                TipStates(ileft) = -2;
            end
            nfcon(1) = ileft;
        end
    end
else
    countp = 1;
    for i = 3: ConnList(Tipele,2)+2
        BndEle = ConnList(Tipele,i);
        if isNF < 0.1
            if BndEle > 0.1 % && Index(BndEle) < 0.1
                if ConnList(BndEle,1) == ConnList(Tipele,1)
                    ConnListnew(Tipele,2) = ConnList(Tipele,2)+1;
                    ConnListnew(Tipele,ConnListnew(Tipele,2)+2) = nAllEle+1;
                    Coord1 = AllEle(Tipele,1:2);
                    isin =  isDotIn(Coord1,AllEle(ConnList(Tipele,i),1:4));
                    if isin > 0.1
                        TipCoord = AllEle(Tipele,3:4);
                        ConnListnew(nAllEle+1:nAllEle+GrowthN,1) = nHF_global + 1;
                        nHF_global = nHF_global +1;
                    else
                        TipCoord = AllEle(Tipele,1:2);
                        ConnListnew(nAllEle+1:nAllEle+GrowthN,1)=nHF_global + 1;
                        nHF_global = nHF_global +1;
                    end
                else
                    if AllEle(BndEle,10) < 2.9
                        TipStates(BndEle) = -2;
                    end
                    nfcon(countp) = BndEle;
                    countp = countp +1 ;
                end
            end
        else
            if BndEle > 0.1 && Index(BndEle) < 0.1
                if ConnList(BndEle,1) == ConnList(Tipele,1)
                    ConnListnew(Tipele,2) = ConnList(Tipele,2)+1;
                    ConnListnew(Tipele,ConnListnew(Tipele,2)+2) = nAllEle+1;
                    Coord1 = AllEle(Tipele,1:2);
                    isin =  isDotIn(Coord1,AllEle(ConnList(Tipele,i),1:4));
                    ConnListnew(nAllEle+1:nAllEle+GrowthN,1) = ConnList(Tipele,1);
                    if isin < 0.1
                        TipCoord = AllEle(Tipele,3:4);
                    else
                        TipCoord = AllEle(Tipele,1:2);
                    end
                    if AllEle(BndEle,10) < 2.9
                        TipStates(BndEle) = -2;
                    end
                    nfcon(countp) = BndEle;
                    countp = countp +1 ;
                end
            end
        end
    end
end

alpha = asind(AllEle(Tipele,5));
beta = alpha + theta;
if beta > 360
    beta = beta - 360;
end

Dirvec = [cosd(beta),sind(beta)];

AllElenew(nAllEle+1:nAllEle+GrowthN,7) = addL;

ConnListnew(nAllEle+1:nAllEle+GrowthN,2) = 2;
ConnListnew(nAllEle+1,3) = Tipele;
ConnListnew(nAllEle+1,4) = nAllEle+2;
if isNF < 0.1
    if ConnList(Tipele,2) > 2.1
        ConnListnew(nAllEle+1,2) = 4;
        ConnListnew(nAllEle+1,5) = nfcon(1);
        ConnListnew(nAllEle+1,6) = nfcon(2);
    end
else
    if ConnList(Tipele,2) > 2.1
        ConnListnew(nAllEle+1,2) = 3;
        ConnListnew(nAllEle+1,5) = nfcon(1);
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

for i = 1 : GrowthN
    x1 = AllElenew(nAllEle+i,1);
    x2 = AllElenew(nAllEle+i,3);
    y1 = AllElenew(nAllEle+i,2);
    y2 = AllElenew(nAllEle+i,4);
    if x1 > x2 ||(x1==x2 && y1 > y2)
        temp = ConnListnew(nAllEle+i,3);
        ConnListnew(nAllEle+i,3) = ConnListnew(nAllEle+i,4);
        ConnListnew(nAllEle+i,4) = temp;
        temp = AllElenew(nAllEle+i,1:2);
        AllElenew(nAllEle+i,1:2) = AllElenew(nAllEle+i,3:4);
        AllElenew(nAllEle+i,3:4) = temp;
    end
end
if beta < -90
    beta = beta + 180;
else
    if beta > 90  && beta <= 270
        beta = beta - 180;
        
    end
    if beta > 270
        beta = beta - 360;
    end
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
Cp = [Cp;0*Cp(1:GrowthN,:)+1e-5*Cp(Index(Tipele),:)];
Xf = [Xf;0*Xf(1:GrowthN,:)+Xf(Index(Tipele),:)];

nAct = nAct + GrowthN;
Index(nAllElenew-GrowthN+1:nAllElenew) = nAct-GrowthN+1:nAct;
IndexInv(nAct-GrowthN+1:nAct) = nAllElenew-GrowthN+1:nAllElenew;

if AllEle(nf,10) < 1.1 && ConnList(Tipele,2) < 2.2
    [iNF,jNF,Point,isMerge,mindis] = CheckCross(Tipele,AllElenew(nAllEle+GrowthN,1:4),nAllEle,AllEle,ConnList,epsDis);
    if isMerge > 0.1
        HasInter  = 1;
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
        
        [nAllElenew,AllElenew,ConnListnew] = BuildNewEle(Point,nAllElenew,AllElenew,ConnListnew,CrossType,iNF,jNF,addL);
        
        le = ConnListnew(nAllEle+GrowthN,3);
        if le < 0.1
            ConnListnew(nAllEle+GrowthN,3) = nAllEle+GrowthN+1;
        else
            ConnListnew(nAllEle+GrowthN,4) = nAllEle+GrowthN+1;
        end
        iparent = ConnList(iNF,1);
        ArrivalNF(iparent)=1;
        %         ShearS = AllEle(iNF,11);
        
        Dshf = Dnew(nAct);
        Dnhf = Dnew(2*nAct);
        Dadd = zeros(GrowthN,2);
        xtip = GrowthN-0.5:-1:0.5;
        xtip = xtip'*addL;
        Li = GrowthN*addL+addL/2;
        Dadd(:,1) = Dshf/sqrt(Li)*sqrt(xtip)*0;
        Dadd(:,2) = Dnhf/sqrt(Li)*sqrt(xtip)*0+Dinit;
        if iNF ~= jNF
            Dnew = [Dnew(1:nAct);Dadd(:,1);
                Dnew(nAct+1:nAct*2);Dadd(:,2)];
            Pnew = [Pnew(1:nAct);zeros(GrowthN,1)+Pinit];
        else
            Dnew = [Dnew(1:nAct);Dadd(:,1);
                Dnew(nAct+1:nAct*2);Dadd(:,2)];
            Pnew = [Pnew(1:nAct);zeros(GrowthN,1)+Pinit];
        end
        hbanking = [hbanking;0*hbanking(1:GrowthN,:)];
        hslurry = [hslurry; 0*hslurry(1:GrowthN,:)+Mat.h];
        Cp = [Cp;0*Cp(1:GrowthN,:)+1e-5*Cp(Index(Tipele),:)];
        Xf = [Xf;0*Xf(1:GrowthN,:)+Xf(Index(Tipele),:)];
        TipStates(nAllEle+GrowthN*2)= 999+state0 ;
        TipStates(nAllEle+GrowthN) = 0;
        nAct = nAct + GrowthN;
        Index(nAllEle+GrowthN+1:nAllEle+GrowthN*2) = nAct-GrowthN+1:nAct;
        IndexInv(nAct-GrowthN+1:nAct) = nAllEle+GrowthN+1:nAllEle+GrowthN*2;
    end
end
AllEle_global(1:nAllElenew,:) = AllElenew;
nAllEle_global = nAllElenew;
ConnList_global(1:nAllElenew,:)  = ConnListnew;
for i = 1 : nAct
    CpF_global(IndexInv(i),:) = Cp(i,:);
    XfF_global(IndexInv(i),:) = Xf(i,:);
    Hbanking_global(IndexInv(i),:) = hbanking(i,:);
    Hslurry_global(IndexInv(i),:) = hslurry(i,:);
    DD_global(IndexInv(i)) = Dnew(i);
    DD_global(IndexInv(i)+MaxEle) = Dnew(i+nAct);
    PresF_global(IndexInv(i)) = Pnew(i);
end
end