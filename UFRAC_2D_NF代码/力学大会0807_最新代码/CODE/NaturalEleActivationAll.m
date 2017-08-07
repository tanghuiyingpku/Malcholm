function    isGrow = NaturalEleActivationAll()
%Normal Stresses
global TipStates Mat AllEle_global  ConnList_global PresF_global DD_global
global nTip MaxEle Tipcoordinate GrowNumber;
global IndexInv nAct Index  isMechActive_global XfF_global  nAllEle_global;
DD_pre = zeros(nAct*2,1);
for i = 1 : nAct
    DD_pre(i) = DD_global(IndexInv(i));
    DD_pre(i+nAct) = DD_global(MaxEle+IndexInv(i));
end
nAct0 = nAct;
Actnum = zeros(50,1);
ActTip = zeros(50,1);
countN = 10;
countT = 0;
AlreadyAct = zeros(nAllEle_global,1);
isGrow = 0;

while countN > 0.1
countN = 0
for i = 1 : nAct0
    if AllEle_global(IndexInv(i),10) - 2 >= 0
        % Nf
        Index_Frac = ConnList_global(IndexInv(i),1);
        base = IndexInv(i);
        Xfinit = XfF_global(IndexInv(i),:);
        nconn = ConnList_global(IndexInv(i),2);
        % 如有其它天然裂缝与当前裂缝相交，出现新的尖端单元
        if nconn > 1.1
            if TipStates(IndexInv(i)) > 0.1 
                for ib = 1 : ConnList_global(base,2)
                    nb = ConnList_global(base,2+ib);
                    if nb < 0.1
                        continue;
                    end
                    jump =0;
                    front = 0;
                    if  Index(nb) < 0.1 && AlreadyAct(nb) < 0.1
                        countN = countN + 1;
                        Actnum(countN) = nb;
                        if ConnList_global(nb,1) == Index_Frac
                            front = nb;
                        else
                            countT = countT+ 1;
                            ActTip(countT) = nb;
                        end
                        AlreadyAct(nb) = 1;
                    end
                    if front < 0.1
                        continue;
                    end
                    ic = front;
                    past = base;
                    count = 0;
                    while  count <  GrowNumber
                        if ic ~= past && ic ~= front
%                             countN = countN + 1;
%                             Actnum(countN) = ic;
                           
                        end
                        if min(ConnList_global(ic,3:6)) < -0.1
                            countT = countT+ 1;
                            ActTip(countT) = ic;
                            jump = 1;
                            break;
                        end
                        for ie = 1 : ConnList_global(ic,2)
                            bc = ConnList_global(ic,2+ie);
                            if bc > 0.1 && bc ~= past
                                if ConnList_global(bc,1) == Index_Frac
                                    %Activate current
                                    if Index(bc) < 0.1 && AlreadyAct(bc) < 0.1
                                        count = count + 1;
                                        countN = countN + 1;
                                        Actnum(countN) = bc;
                                        AlreadyAct(bc) = 1;
                                    end
                                    past = ic;
                                    ic = bc;
                                else
                                    %cross with other natural fractures
                                    if Index(bc) < 0.1 && AlreadyAct(bc) < 0.1
                                        if bc == 37
                                            f = 1;
                                        end
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
                    if jump < 0.1
                        countT = countT+ 1;
                        ActTip(countT) = ic;
                    end
                    TipStates(IndexInv(i)) = 0;
                    nTip = nTip-1;
                end
            end
        end
    end
end
for i = 1:countN
    nb = Actnum(i);
    AllEle_global(nb,10) = 3;
    PresF_global(nb) = Mat.Pp;
    XfF_global(nb,:) = Xfinit;
    Index(nb) = nAct+1;
    IndexInv(nAct+1) = nb;
    nAct = nAct + 1;
    isMechActive_global(nAct) = -2;
end
for i = 1 : countT
    bc = ActTip(i);
    TipStates(bc) = nTip+1;
    Tipcoordinate(bc,1:2) = FindTipCoord(bc);
    nTip = nTip + 1;
end
if nAct > nAct0
    isGrow = 1;
end
end
clear DD_pre Actnum ActTip AlreadyAct;
end
