function  DynamicGrowthNF(Pinit,Dinit,nf)
global TipStates MaxEle DD_global  CpF_global XfF_global Tipcoordinate ArrivalNF;
global IndexInv nAct isMechActive_global Index  ConnList_global AllEle_global PresF_global nTip;

state0 = TipStates(nf);
TipStates(nf) = 0;
%nTip = nTip - 1;
%nf 裂缝编号
%index 1,2 代表是哪一个tip
%增加两个element
%新增网格编号
for i = 1 : 1
    ncon = ConnList_global(nf,2);
    for j = 1 : ncon
        nb = ConnList_global(nf,j+2);
        if nb < 0.1
            continue;
        end
        if TipStates(nb) > 0.1 && ConnList_global(nb,2) < 2.1 && min(ConnList_global(nb,3:4)) > 0.1
            % Two tip emerge
            TipStates(nb) = 0;
            nTip = nTip -2;
             continue;
        end
        if Index(nb) < 0.1%AllEle_global(nb,10) > 1.1
            if 1;%abs(AllEle_global(nb,10) -3) > 1e-3
                ic = nb;
                %nTip = nTip + 1;
                if ConnList_global(ic,1) == ConnList_global(nf,1) && AllEle_global(ic,10) < 2.9
                    nAct = nAct + 1;
                    IndexInv(nAct) = ic;
                    Index(ic) = nAct;
                    isMechActive_global(nAct) = -2;
                    PresF_global(ic) = Pinit/2;
                    DD_global(ic) = DD_global(nf)*1e-5;
                    if abs(DD_global(nf)) < 1e-16 && ConnList_global(ic,2) < 2.2
                        DD_global(ic) = DD_global(ConnList_global(nf,3-j+2))*1e-5;
                    end
                    DD_global(ic+MaxEle) = Dinit;
                    CpF_global(ic,:) = CpF_global(nf,:);
                    XfF_global(ic,:) = XfF_global(nf,:);
                    AllEle_global(nb,10) = 3;
                    TipStates(ic) = state0;
                    point1 = AllEle_global(ic,1:2);
                    isIn = isDotIn(point1,AllEle_global(nf,1:4));
                    if isIn < 0.1
                        Tipcoordinate(ic,:) = point1;
                    else
                        Tipcoordinate(ic,:) = AllEle_global(ic,3:4);
                    end
                end
                if ConnList_global(ic,1) ~= ConnList_global(nf,1) % Add new tip at another NF
                    if AllEle_global(nb,10) > 1.1
                        ArrivalNF(ConnList_global(ic,1) )=1;
                    end
                    nTip = nTip + 1;
                    nAct = nAct + 1;
                    IndexInv(nAct) = ic;
                    Index(ic) = nAct;
                    isMechActive_global(nAct) = -2;
                    PresF_global(ic) = Pinit/2;
                    DD_global(ic) = DD_global(nf)*1e-5;
                     if abs(DD_global(nf)) < 1e-16
                        DD_global(ic) = DD_global(ConnList_global(nf,3-j+2))*1e-5;
                    end
                    DD_global(ic+MaxEle) = Dinit;
                    CpF_global(ic,:) = CpF_global(nf,:);
                    XfF_global(ic,:) = XfF_global(nf,:);
                    if AllEle_global(nb,10) > 1.1
                        AllEle_global(nb,10) = 3;
                    end
                    TipStates(ic) = nTip;
                    point1 = AllEle_global(ic,1:2);
                    isIn = isDotIn(point1,AllEle_global(nf,1:4));
                    if isIn < 0.1
                        Tipcoordinate(ic,:) = point1;
                    else
                        Tipcoordinate(ic,:) = AllEle_global(ic,3:4);
                    end
                end
            end
        end
        
    end
end
end