function [nEle2,Ele2,ConnList2] = BuildNewEle_backup(Point,nAllEle,AllEle,ConnList,CrossType,iNF,jNF,Lseq)
%Type 1 : 裂缝相交在某个单元中间
%Type 2: 裂缝相交在某个单元的端点
global  Index DD_global MaxEle PresF_global Tipcoordinate;
global nAct
if nAct > 120
    f = 1;
end
GrowthHF = length(Lseq);
maxN = 3;
minN = 2;
addL = Lseq(1);
Ds=0;
Dn=0;
Pinit =0;
if Index(iNF) > 0.1
    Ds =  DD_global(iNF);
    Dn =  DD_global(iNF+ MaxEle);
    Pinit = PresF_global(iNF);
end
if Index(jNF) > 0.1
    Ds =  DD_global(jNF);
    Dn =  DD_global(jNF+ MaxEle);
    Pinit = PresF_global(jNF);
end

if CrossType < 1.9  ||   jNF < 0.1 %At the tip of NF
    right = ConnList(iNF,4);
    NfEle = AllEle(iNF,1:4);
    dis1 = CalculateDis(NfEle(1:2),Point);
    dis2 = CalculateDis(NfEle(3:4),Point);
    if jNF < 0.1
        nL = maxN;
        nR = minN;
    else
        nL = ceil(dis1/addL);
        if nL < minN
            nL = minN;
        end
        if nL > maxN
            nL  = maxN;
        end
        nR = ceil(dis2/addL);
        if nR < minN
            nR = minN;
        end
        if nR > maxN
            nR = maxN;
        end
    end
    addN = nL + nR;
    %加上两个连接的压裂裂缝，减去一个现有的天然裂缝编号
    nEle2 = nAllEle+addN + GrowthHF -1;
    %设置空间
    Ele2 = [AllEle;zeros(nEle2-nAllEle,11)];
    [~,nc] = size(ConnList);
    ConnList2 = [ConnList;zeros(nEle2-nAllEle,nc)];
    
    %修改天然裂缝的部分
    Ele2(nAllEle+GrowthHF+1:nEle2,5) = AllEle(iNF,5);
    Ele2(nAllEle+GrowthHF+1:nEle2,6) = AllEle(iNF,6);
    Ele2(nAllEle+GrowthHF+1:nEle2,10) = AllEle(iNF,10);
    Ele2(nAllEle+GrowthHF+1:nEle2,11) = AllEle(iNF,11);
    Ele2(iNF,3:4) =Ele2(iNF,1:2)+(Point-Ele2(iNF,1:2))/nL;
    DD_global(nAllEle+GrowthHF+1:nEle2) = Ds;
    DD_global(nAllEle+GrowthHF+1+MaxEle:nEle2+MaxEle) = Dn;
    PresF_global(nAllEle+GrowthHF+1:nEle2)=Pinit;
    
    ConnList2(iNF,4) = nAllEle+GrowthHF+1;
    ConnList2(nAllEle+GrowthHF+1,2) = 2;
    ConnList2(nAllEle+GrowthHF+1,3) = iNF;
    if nL + nR > 2.1
        ConnList2(nAllEle+GrowthHF+1,4) = nAllEle+ GrowthHF + 2;
    end
    for i = 1 : nL-1
        Ele2(nAllEle+GrowthHF+i,1:2) = Ele2(iNF,1:2)+(Point-Ele2(iNF,1:2))/nL*i;
        Ele2(nAllEle+GrowthHF+i,3:4) = Ele2(iNF,1:2)+(Point-Ele2(iNF,1:2))/nL*(i+1);
        if i >= 2
            ConnList2(nAllEle+GrowthHF+i,2) = 2;
            ConnList2(nAllEle+GrowthHF+i,3) =nAllEle+ GrowthHF+i-1;
            ConnList2(nAllEle+GrowthHF+i,4) =nAllEle+ GrowthHF+i+1;
        end
    end
    if nL < 1.1
        i = 0;
    end
    if nL > 1.1
        nf1 = nAllEle+GrowthHF+i;
        nf2 = nAllEle+GrowthHF+i+1;
        Ele2(nAllEle+GrowthHF+i,3:4) = Point;
    else
        nf1 = iNF;
        nf2 = nAllEle+GrowthHF+1;
    end
    
    %还是保持为2的类型，因为不确定哪个方位先打开！！！！！
    %Ele2(nAllEle+GrowthHF+i,10) = 3;
    %Ele2(nAllEle+GrowthHF+i,10) = 3;
    %母裂缝都是一样的
    ConnList2(nAllEle+GrowthHF+1:nEle2,1) = ConnList2(iNF,1);
    %裂缝类型改为Open的天然裂缝
    if nL > 1.1
        ConnList2(nAllEle+GrowthHF+i,2) = 3;
        if i > 1
            ConnList2(nAllEle+GrowthHF+i,3) = nAllEle+GrowthHF+i-1;
        else
            ConnList2(nAllEle+GrowthHF+i,3) = iNF;
        end
        ConnList2(nAllEle+GrowthHF+i,4) = nAllEle+GrowthHF+i+1;
        ConnList2(nAllEle+GrowthHF+i,5) = nAllEle+GrowthHF;
    else
        ConnList2(nAllEle+GrowthHF+1,2) = 3;
        ConnList2(nAllEle+GrowthHF+1,3) = iNF;
        ConnList2(nAllEle+GrowthHF+1,4) = nAllEle+GrowthHF+2;
        ConnList2(nAllEle+GrowthHF+1,5) = nAllEle+GrowthHF;
    end
    %与压裂裂缝相邻的左右两个网格
    %     TipStates(nAllEle+GrowthHF+i) = nTip+1;
    %     TipStates(nAllEle+GrowthHF+i+1) = nTip+2;
    %     nTip = nTip + 2;
    %改变单元类型
    %还是保持为2的类型，因为不确定哪个方位先打开！！！！！
    %Ele2(nAllEle+GrowthHF+i+1,10) = 3;
    if nL > 1.1
        ConnList2(nAllEle+GrowthHF+i+1,2) =3;
        ConnList2(nAllEle+GrowthHF+i+1,3) = nAllEle+GrowthHF+i;
        if nAllEle+GrowthHF+i+2 > nEle2
            ConnList2(nAllEle+GrowthHF+i+1,4) = right;
        else
            ConnList2(nAllEle+GrowthHF+i+1,4) = nAllEle+GrowthHF+i+2;
        end
        ConnList2(nAllEle+GrowthHF+i+1,5) = nAllEle+GrowthHF;
    end
    ii = 0;
    if nL < 1.1
        Ele2(nAllEle+GrowthHF+ii+1,1:2) = Point;
    end
    if nL > 1.1
    for ii = nAllEle+GrowthHF+i+1: nEle2
        if ii > nAllEle+GrowthHF+i+1 && ii < nEle2
            ConnList2(ii,2)= 2;
            ConnList2(ii,3) = ii-1;
            ConnList2(ii,4) = ii+1;
        end
        Ele2(ii,1:2) = Ele2(ii-1,3:4);
        Ele2(ii,3:4) = Point+(NfEle(3:4)-Point)/nR*(ii-nAllEle-GrowthHF-i);
    end
    end
   
    Ele2(nEle2,3:4) = NfEle(3:4);
    
    if nAllEle+GrowthHF+i+2 <= nEle2
        ConnList2(nEle2,2) = 2;
    end
    if nL + nR < 2.1
        ConnList2(nEle2,3) = iNF;
    else
        ConnList2(nEle2,3) = nEle2-1;
    end
    ConnList2(nEle2,4)=right;
    if right > 0.1
        ConnList2(right,3) = nEle2;
    end
    
    Ele2(iNF,8) = 0.5*(Ele2(iNF,1)+Ele2(iNF,3));
    Ele2(iNF,9) = 0.5*(Ele2(iNF,2)+Ele2(iNF,4));
    Ele2(iNF,7) = CalculateDis([Ele2(iNF,1) Ele2(iNF,2)],[Ele2(iNF,3),Ele2(iNF,4)]);
    %压裂裂缝上的新单元
    Ele2(nAllEle+1:nAllEle+GrowthHF,10) = 1;
    
    sint = AllEle(nAllEle,5);
    cost = AllEle(nAllEle,6);
    Ele2(nAllEle+1:nAllEle+GrowthHF,5) =sint;
    Ele2(nAllEle+1:nAllEle+GrowthHF,6) = cost;
    Tipele = nAllEle;
    iright = ConnList(Tipele,4);
    if iright < 0.1
        TipCoord = AllEle(Tipele,3:4);
    else
        TipCoord = AllEle(Tipele,1:2);
    end
    TipCoord = Tipcoordinate(Tipele,:);
    nH = GrowthHF;
    HFele = Point - TipCoord;
    vec2 = HFele/sqrt(HFele(1)^2+HFele(2)^2);
  %  Lseq = fliplr(Lseq')';
    Lseq = [0;Lseq];
    for i = 1 : nH
        Ele2(nAllEle+i,1:2) = TipCoord + vec2*sum(Lseq(1:i));
        Ele2(nAllEle+i,3:4) = TipCoord + vec2*sum(Lseq(1:i+1));
        ConnList2(nAllEle+i,2) = 2;
        ConnList2(nAllEle+i,3) = nAllEle+(i-1);
        ConnList2(nAllEle+i,4) = nAllEle+(i+1);
    end
    Ele2(nAllEle+i,3:4) = Point;
    ConnList2(nAllEle+i,2) = 3;
    ConnList2(nAllEle+i,3) = nAllEle+i-1;
    ConnList2(nAllEle+i,4) = nf1;
    ConnList2(nAllEle+i,5) = nf2;
    ConnList2(nAllEle+1:nAllEle+GrowthHF,1) = ConnList(nAllEle,1);
    
    
    vec = [cost,sint];
    if dot(vec,vec2) < 0
        %与天然裂缝相交的单位就不交换了
        for ii = 1 : GrowthHF
            temp = ConnList2(nAllEle+ii,3);
            ConnList2(nAllEle+ii,3) = ConnList2(nAllEle+ii,4);
            ConnList2(nAllEle+ii,4) = temp;
            temp2 = Ele2(nAllEle+ii,1:2);
            Ele2(nAllEle+ii,1:2) = Ele2(nAllEle+ii,3:4);
            Ele2(nAllEle+ii,3:4) = temp2;
        end
        %
    end
    for ii = 1 : nEle2-nAllEle
        Ele2(nAllEle+ii,8) = 0.5*(Ele2(nAllEle+ii,1)+Ele2(nAllEle+ii,3));
        Ele2(nAllEle+ii,9) = 0.5*(Ele2(nAllEle+ii,2)+Ele2(nAllEle+ii,4));
        Ele2(nAllEle+ii,7) = CalculateDis([Ele2(nAllEle+ii,1) Ele2(nAllEle+ii,2)],[Ele2(nAllEle+ii,3),Ele2(nAllEle+ii,4)]);
    end
else
    %如果直接相交于某一点，该点为iNF和jNF间的交点
    %加密iNF和jNF
    Li = CalculateDis(AllEle(iNF,1:2),Point);
    Lj = CalculateDis(AllEle(jNF,3:4),Point);
    %         Li = AllEle(iNF,7);
    %         Lj = AllEle(jNF,7);
    nL = ceil(Li/addL);
    nR = ceil(Lj/addL);
    if nL < minN
        nL = minN;
    end
    if nL > maxN
        nL  = maxN;
    end
    if nR < minN
        nR = minN;
    end
    if nR > maxN
        nR = maxN;
    end
    addN = nL + nR;
    %加上两个连接的压裂裂缝，减去2个现有的天然裂缝编号!
    nEle2 = nAllEle+addN + GrowthHF -2;
    DD_global(nAllEle+GrowthHF+1:nEle2) = Ds;
    DD_global(nAllEle+GrowthHF+1+MaxEle:nEle2+MaxEle) = Dn;
    PresF_global(nAllEle+GrowthHF+1:nEle2)=Pinit;
    
    %设置空间
    Ele2 = [AllEle;zeros(nEle2-nAllEle,11)];
    [~,nc] = size(ConnList);
    ConnList2 = [ConnList;zeros(nEle2-nAllEle,nc)];
    %左侧天然裂缝
    Ele2(iNF,3:4) =Ele2(iNF,1:2)+(Point-Ele2(iNF,1:2))/nL;
    if nL > 1.1
        Ele2(nAllEle+GrowthHF+1:nAllEle+GrowthHF+nL-1,5) = AllEle(iNF,5);
        Ele2(nAllEle+GrowthHF+1:nAllEle+GrowthHF+nL-1,6) = AllEle(iNF,6);
        Ele2(nAllEle+GrowthHF+1:nAllEle+GrowthHF+nL-1,10) = AllEle(iNF,10);
        Ele2(nAllEle+GrowthHF+1:nAllEle+GrowthHF+nL-1,11) = AllEle(iNF,11);
    end
    if nL > 1.1
        Ele2(nAllEle+GrowthHF+nL:nEle2,10) = AllEle(jNF,10);
        Ele2(nAllEle+GrowthHF+nL:nEle2,5) = AllEle(jNF,5);
        Ele2(nAllEle+GrowthHF+nL:nEle2,6) =AllEle(jNF,6);
        Ele2(iNF,3:4) =Ele2(iNF,1:2)+(Point-Ele2(iNF,1:2))/nL;
        ConnList2(iNF,4) = nAllEle+GrowthHF+1;
        ConnList2(nAllEle+GrowthHF+1,2) = 2;
        ConnList2(nAllEle+GrowthHF+1,3) = iNF;
        ConnList2(nAllEle+GrowthHF+1,4) = nAllEle+ GrowthHF + 2;
    else
        ConnList2(iNF,2) = 3;
        ConnList2(iNF,5) = nAllEle+GrowthHF;
    end
    for i = 1 : nL-1
        Ele2(nAllEle+GrowthHF+i,1:2) = Ele2(iNF,1:2)+(Point-Ele2(iNF,1:2))/nL*i;
        Ele2(nAllEle+GrowthHF+i,3:4) = Ele2(iNF,1:2)+(Point-Ele2(iNF,1:2))/nL*(i+1);
        if i >= 2
            ConnList2(nAllEle+GrowthHF+i,2) = 2;
            ConnList2(nAllEle+GrowthHF+i,3) =nAllEle+ GrowthHF+i-1;
            ConnList2(nAllEle+GrowthHF+i,4) =nAllEle+ GrowthHF+i+1;
        end
    end
    if nL < 1.1
        i = 0;
    end
    if nL < 1.1
        nf1 = iNF;
    else
        nf1 = nAllEle+GrowthHF+i;
    end
    nf2 = jNF;
    %压裂裂缝相交的小段开裂
    %还是保持为2的类型，因为不确定哪个方位先打开！！！！！
    %Ele2(nAllEle+GrowthHF+i,10) = 3;
    %母裂缝都是一样的
    if nL > 1.1
        Ele2(nAllEle+GrowthHF+i,3:4) = Point;
        ConnList2(nAllEle+GrowthHF+1:nEle2,1) = ConnList2(iNF,1);
    end
    
    %新增连接
    ConnList2(nAllEle+GrowthHF+i,2) = 3;
    ConnList2(nAllEle+GrowthHF+i,3) = nAllEle+GrowthHF+i-1;
    if nL > 1.1
        if i < 1.1
            ConnList2(nAllEle+GrowthHF+i,3) = iNF;
        else
            ConnList2(nAllEle+GrowthHF+i,3) = nAllEle+GrowthHF+i-1;
        end
        ConnList2(nAllEle+GrowthHF+i,4) = jNF;
        ConnList2(nAllEle+GrowthHF+i,5) = nAllEle+GrowthHF;
    end
    
    %与压裂裂缝相邻的左右两个网格
    %TipStates(nAllEle+GrowthHF+i) = nTip+1;
    %TipStates(jNF) = nTip+2;
    %nTip = nTip + 2;
    RTip = AllEle(jNF,3:4);
    %单元类型的修改
    right = ConnList2(jNF,4);
    %还是保持为2的类型，因为不确定哪个方位先打开！！！！！
    %Ele2(jNF,10) = 3;
    Ele2(jNF,1:2) = Point;
    Ele2(jNF,3:4) = Point+(RTip-Point)/nR;
    %右侧天然裂缝
    %新增连接
    ConnList2(jNF,2) =3;
    if nL > 1.1
        ConnList2(jNF,3) = nAllEle+GrowthHF+i;
    end
    if nR > 1.1
        ConnList2(jNF,4) = nAllEle+GrowthHF+i+1;
    end
    ConnList2(jNF,5) = nAllEle+GrowthHF;
    if nR > 1.1
        ii= nAllEle+GrowthHF+i+1;
        ConnList2(ii,2)= 2;
        ConnList2(ii,3) = jNF;
        ConnList2(ii,4) = ii+1;
        Ele2(ii,1:2) = Ele2(jNF,3:4);
        Ele2(ii,3:4) = Point+(RTip-Point)/nR*2;
        if ii >= nEle2
            ConnList2(nEle2,3) = jNF;
        else
            for ii = nAllEle+GrowthHF+i+2: nEle2
                if ii > nAllEle+GrowthHF+i+1 && ii < nEle2
                    ConnList2(ii,2)= 2;
                    ConnList2(ii,3) = ii-1;
                    ConnList2(ii,4) = ii+1;
                    Ele2(ii,1:2) = Ele2(ii-1,3:4);
                    Ele2(ii,3:4) = Point+(RTip-Point)/nR*(ii-nAllEle-GrowthHF-i+1);
                end
            end
            Ele2(ii,1:2) = Ele2(ii-1,3:4);
            Ele2(ii,3:4) = RTip;
            ConnList2(nEle2,3) = nEle2-1;
        end
        ConnList2(nEle2,2) = 2;
        ConnList2(nEle2,4) = right;
        if right > 0.1
            ConnList2(right,3) = nEle2;
        end
    end
    
    Ele2(iNF,8) = 0.5*(Ele2(iNF,1)+Ele2(iNF,3));
    Ele2(iNF,9) = 0.5*(Ele2(iNF,2)+Ele2(iNF,4));
    Ele2(iNF,7) = CalculateDis([Ele2(iNF,1) Ele2(iNF,2)],[Ele2(iNF,3),Ele2(iNF,4)]);
    
    Ele2(jNF,8) = 0.5*(Ele2(jNF,1)+Ele2(jNF,3));
    Ele2(jNF,9) = 0.5*(Ele2(jNF,2)+Ele2(jNF,4));
    Ele2(jNF,7) = CalculateDis([Ele2(jNF,1) Ele2(jNF,2)],[Ele2(jNF,3),Ele2(jNF,4)]);
    
    %压裂裂缝上的新单元
    Ele2(nAllEle+1:nAllEle+GrowthHF,10) = 1;
    
    sint = AllEle(nAllEle,5);
    cost = AllEle(nAllEle,6);
    Ele2(nAllEle+1:nAllEle+GrowthHF,5) =sint;
    Ele2(nAllEle+1:nAllEle+GrowthHF,6) = cost;
    Tipele = nAllEle;
    iright = ConnList(Tipele,4);
    if iright < 0.1
        TipCoord = AllEle(Tipele,3:4);
    else
        TipCoord = AllEle(Tipele,1:2);
    end
    TipCoord = Tipcoordinate(Tipele,:);
    nH = GrowthHF;
    HFele = Point - TipCoord;
    vec2 = HFele/sqrt(HFele(1)^2+HFele(2)^2);
  %  Lseq = fliplr(Lseq')';
    Lseq = [0;Lseq];
    for i = 1 : nH
        Ele2(nAllEle+i,1:2) = TipCoord + vec2*sum(Lseq(1:i));
        Ele2(nAllEle+i,3:4) = TipCoord + vec2*sum(Lseq(1:i+1));
        ConnList2(nAllEle+i,2) = 2;
        ConnList2(nAllEle+i,3) = nAllEle+(i-1);
        ConnList2(nAllEle+i,4) = nAllEle+(i+1);
    end
    Ele2(nAllEle+i,3:4) = Point;
    ConnList2(nAllEle+i,2) = 3;
    ConnList2(nAllEle+i,3) = nAllEle+i-1;
    ConnList2(nAllEle+i,4) = nf1;
    ConnList2(nAllEle+i,5) = nf2;
    ConnList2(nAllEle+1:nAllEle+GrowthHF,1) = ConnList(nAllEle,1);
    
    vec = [cost,sint];
    if dot(vec,vec2) < 0
        %与天然裂缝相交的单位就不交换了
        for ii = 1 : GrowthHF
            temp = ConnList2(nAllEle+ii,3);
            ConnList2(nAllEle+ii,3) = ConnList2(nAllEle+ii,4);
            ConnList2(nAllEle+ii,4) = temp;
            temp2 = Ele2(nAllEle+ii,1:2);
            Ele2(nAllEle+ii,1:2) = Ele2(nAllEle+ii,3:4);
            Ele2(nAllEle+ii,3:4) = temp2;
        end
        %
    end
    for ii = 1 : nEle2-nAllEle
        Ele2(nAllEle+ii,8) = 0.5*(Ele2(nAllEle+ii,1)+Ele2(nAllEle+ii,3));
        Ele2(nAllEle+ii,9) = 0.5*(Ele2(nAllEle+ii,2)+Ele2(nAllEle+ii,4));
        Ele2(nAllEle+ii,7) = CalculateDis([Ele2(nAllEle+ii,1) Ele2(nAllEle+ii,2)],[Ele2(nAllEle+ii,3),Ele2(nAllEle+ii,4)]);
    end
end
end