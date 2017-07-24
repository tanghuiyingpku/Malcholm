function KinkEleActivation(dt,CurT)
%Normal Stresses
global TipStates Mat AllEle_global nAllEle_global ConnList_global PresF_global DD_global
global nTip MaxEle Tipcoordinate;
global IndexInv nAct Index CpF_global XfF_global ;
% Use Virtual Element Method
nEle = nAct + 1; % Active one element each time
DD = zeros(nEle*2,1);
Pres = zeros(nEle,1);
AllEle = zeros(nEle,11);
ConnList = zeros(nEle,8);
% Set values for the existing elements
for i = 1 : nAct
    DD(i) = DD_global(IndexInv(i));
    DD(i+nEle) = DD_global(MaxEle+IndexInv(i));
    Pres(i) = PresF_global(IndexInv(i));
    AllEle(i,:)  = AllEle_global(IndexInv(i),:);
    ConnList(i,:) = ConnList_global(IndexInv(i),:);
end

for i = 1 : nAllEle_global
    if abs(TipStates(i)+3) < 1e-3
        %The elements been crossed by HF
        nconn = ConnList_global(i,2);
        for ii = 1 : nconn
            bnd = ConnList_global(i,2+ii);
            if bnd > -0.1
                if abs(AllEle_global(bnd,10) -1) < 1e-3
					Pinit = PresF_global(bnd);
                    Pres(nEle) = Pinit;
					AllEle(nEle,:) = AllEle_global(i,:);
					ConnList(nEle,:) = ConnList_global(i,:);
                    DD_test= CalcDD_Cross(dt,CurT,nEle,DD,Pres,AllEle,ConnList);
                    Ds = DD_test(nEle);
					%DD_test = CalcDD_test(Mat,nEle,Pres,AllEle);
                    % Get KIC at this nf element
					Dn = -DD_test(nEle*2);
                    d = AllEle(nEle,7);
                    K1 = -Mat.G*1e3/4/(1-Mat.miu)*(2*pi/d)^0.5*Dn;
                    K2 = abs( -Mat.G*1e3/4/(1-Mat.miu)*(2*pi/d)^0.5*Ds);
                    Gi = (K1/Mat.K1NF)^2 + (K2/Mat.K2NF)^2;
                    if Gi > 1 %&& K1 > 1e-6
                        point1 = AllEle_global(i,1:2);
                        isIn = isDotIn(point1,AllEle_global(bnd,1:4));
                        if isIn < 0.1
                            Tipcoordinate(i,:) = point1;
                        else
                            Tipcoordinate(i,:) = AllEle_global(i,3:4);
                        end
                        nAct = nAct+1;
                        Index(i) = nAct;
                        IndexInv(nAct) = i;
                        %Activate i
                        PresF_global(IndexInv(nAct)) = Pinit;
                        DD_global(IndexInv(nAct)) = 0;
                        CpF_global(IndexInv(nAct),:) = CpF_global(bnd,:);
                        XfF_global(IndexInv(nAct),:) = XfF_global(bnd,:);
                        DD_global(IndexInv(nAct)+MaxEle) = 0;
                        AllEle_global(i,10) = 5;
                        nTip = nTip + 1;
                        TipStates(i) = nTip;
                        break;
                    end
                end
                %end
            end
        end
    end
end
clear DD_test Pres DD ConnList AllEle;
end
