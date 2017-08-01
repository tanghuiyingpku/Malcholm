function  Precondition_NF(CurT, Elei0)
global Mat  ActEle nAct 
global PresF_global  IndexInv CpF_global XfF_global ;
global Index   Tau_global TipStates nTip;

Pres = Mat.Pp;
Elej = Elei0(2);% ConnList_global(Elei,3);
Elei = Elei0(1);

nTip = nTip+ 1;
CpF_global(Elei,:) = 0;
XfF_global(Elei,1) = 1;
ActEle(Elei) = 1;
PresF_global(Elei) = Pres;
Tau_global(Elei) = CurT;
IndexInv(nAct+1) = Elei;
Index(Elei) = nAct+1;
TipStates(Elei) = nTip;

nTip = nTip+ 1;

ActEle(Elej) = 1;
XfF_global(Elej,1) = 1;
PresF_global(Elej) =  Pres;
Tau_global(Elej) = CurT;
IndexInv(nAct+2) = Elej;
Index(Elej) = nAct+2;
CpF_global(Elej,:) = 0;
nAct = nAct + 2;
TipStates(Elej) = nTip;
end