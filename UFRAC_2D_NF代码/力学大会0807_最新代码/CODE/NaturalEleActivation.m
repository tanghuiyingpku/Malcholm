function    NaturalEleActivation()
%Normal Stresses
global TipStates Mat AllEle_global  ConnList_global PresF_global DD_global
global nTip MaxEle Tipcoordinate;
global IndexInv nAct Index  ArrivalNF isMechActive_global XfF_global Fractures;
DD_pre = zeros(nAct*2,1);
for i = 1 : nAct
    DD_pre(i) = DD_global(IndexInv(i));
    DD_pre(i+nAct) = DD_global(MaxEle+IndexInv(i));
end
nAct0 = nAct;
for i = 1 : nAct0
    if AllEle_global(IndexInv(i),10) - 2 >= 0
        % Nf
        Index_Frac = ConnList_global(IndexInv(i),1);
        nconn = ConnList_global(IndexInv(i),2);
        % 如有其它天然裂缝与当前裂缝相交，出现新的尖端单元
        if nconn > 2.1
            if abs(DD_pre(i)) > 1e-6 || abs(DD_pre(i+nAct0)) > 1e-6 || abs(PresF_global(IndexInv(i)) - Mat.Pp)>Mat.Pp*0.15
                for ib = 1 : ConnList_global(IndexInv(i),2)
                    nb = ConnList_global(IndexInv(i),2+ib);
                    if nb < 0.1
                        continue;
                    end
                    Index_Frac2 = ConnList_global(nb,1);
                    if  Index_Frac2 ~= Index_Frac
                        if Index(nb)< 0.1
                            %这条裂缝还没有流体进入，则激活当前裂缝
                            ArrivalNF(Index_Frac2) = 1;
                            for ie2 = Fractures{Index_Frac2}.EleFirst:Fractures{Index_Frac2}.EleLast
                                Xfinit = XfF_global(IndexInv(i),:);
                                AllEle_global(ie2,10) = 3;
                                PresF_global(ie2) = Mat.Pp;
                                XfF_global(ie2,:) = Xfinit;
                                Index(ie2) = nAct+1;
                                IndexInv(nAct+1) =ie2 ;
                                nAct = nAct + 1;
                                isMechActive_global(nAct) = -2;
                                if min(ConnList_global(ie2,3:5)) < -0.1
                                    TipStates(ie2) = nTip+1;
                                    Tipcoordinate(ie2,1:2) = FindTipCoord(ie2);
                                    nTip = nTip + 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
clear DD_pre;
end
