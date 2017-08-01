function DynamicGrowth_2(iscross, Pinit,Dinit,nf,theta)
global ArrivalNF;
global stepL Mat;
global GrowthN Tipcoordinate;
global HasInter MaxEle;
global TipStates Hslurry_global Hbanking_global CpF_global XfF_global;
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

epsDis = stepL*1.2;
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
            if AllEle(BndEle,10) > 1.1
                nfcon(countp) = BndEle;
                countp = countp +1 ;
            end
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
    alpha = 0.7;
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
        AllElenew(nAllEle+i,9) = AllElenew(Tipele,2)*0.5+AllElenew(Tipele,4)*0.5;
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