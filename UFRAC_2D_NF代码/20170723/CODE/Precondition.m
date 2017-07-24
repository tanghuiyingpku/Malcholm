function  Precondition(CurT, Elei,rate)
global Mat  ActEle nAct ConnList_global
global PresF_global  IndexInv CpF_global XfF_global stepL;
global Index nInitial Fluid Tau_global;
L = stepL*nInitial;
xx = -stepL*(nInitial-0.5): stepL : stepL*(nInitial-0.5);
x1 = -stepL*(nInitial): stepL : stepL*(nInitial-1);
x2 = -stepL*(nInitial-1): stepL : stepL*(nInitial);
Pana = zeros(2*nInitial,1);
Dana = zeros(2*nInitial,1);
visco = Fluid.fluid{1}.visco;
Pi = -Mat.Sxx*1.01*1e6;
for i = 1 : nInitial*2
    Pana(i) = (16*visco*rate*2*(Mat.E*1e9)^3/pi/Mat.h^4*(L-abs(xx(i))))^0.25-Mat.Sxx*1e6*1.01;
    Dana(i) = 2*(Pana(i)+Mat.Sxx*1e6)*L*2/Mat.E/1e9;
end
nEle = nInitial*2;
AllEle = zeros(nEle,11);
ConnList = zeros(nEle,8);
AllEle(:,1) = x1;
AllEle(:,3) = x2;
AllEle(:,2) = 0;
AllEle(:,4) = 0;
AllEle(:,5) = 0;
AllEle(:,6) = 1;
AllEle(:,8) = 0.5*(AllEle(:,1)+AllEle(:,3));
AllEle(:,9) = 0.5*(AllEle(:,2)+AllEle(:,4));
AllEle(:,7) = stepL;
AllEle(:,10) = 1;
ConnList(:,1) = 1;
ConnList(:,2) = 2;
for i = 1 : nEle
    ConnList(i,3) = i-1;
    ConnList(i,4) = i+1;
    if i < 1.1
        ConnList(i,3) = -1;
    end
    if i > nEle-0.1
        ConnList(i,4) = -1;
    end
end
% Calculate Initial Pressure Distribution - steady state, no need for
%dt
%[Pres,DD] = Intial_Balance(nEle,AllEle,ConnList,Pana,rate,Dana);
Pres = Pana/1e6;
Elej = Elei+1;% ConnList_global(Elei,3);


CpF_global(Elei,:) = 0;
XfF_global(Elei,1) = 1;
ActEle(Elei) = 1;
PresF_global(Elei) = Pres(nInitial);
Tau_global(Elei) = CurT;
IndexInv(nAct+1) = Elei;
Index(Elei) = nAct+1;

ActEle(Elej) = 1;
XfF_global(Elej,1) = 1;
PresF_global(Elej) =  Pres(nInitial);
Tau_global(Elej) = CurT;
IndexInv(nAct+2) = Elej;
Index(Elej) = nAct+2;
CpF_global(Elej,:) = 0;
nAct = nAct + 2;
for i = 1 : nInitial-1
    for kk =  1 : 2
        num1 = ConnList_global(Elei,kk+2);
        if num1 > 0.1
            if ActEle(num1) < 0.1
                ActEle(num1) = 1;
                XfF_global(num1,1) = 1;
                PresF_global(num1) =  Pres(nInitial-i);
                Tau_global(num1) = CurT;
                CpF_global(num1,:) = 0;
                IndexInv(nAct+1) = num1;
                Index(num1) = nAct+1;
                nAct = nAct + 1;
                i1 = num1;
            end
        end
    end
     for kk =  1 : 2
         num2 = ConnList_global(Elej,kk+2);
         if num2 > 0.1
             if ActEle(num2) < 0.1
                 ActEle(num2) = 1;
                 XfF_global(num2,1) = 1;
                 PresF_global(num2) =  Pres(nInitial-i);
                 CpF_global(num2,:) = 0;
                 Tau_global(num2) = CurT;
                 IndexInv(nAct+1) = num2;
                 Index(num2) = nAct+1;
                 nAct = nAct + 1;
                 j1 = num2;
             end
         end
     end
     Elei = i1;
     Elej = j1;
end
end