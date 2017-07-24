function [DD_new,Pnew,AllElenew,nElenew,ConnListnew] = GrowthPathJudge(EleType,Mat,Pres,DD,dt,nAllEle,AllEle,ConnList,nFracture,Fractures,Sxx,Syy,Sxy,nwell,well)
global KIChf KICnf
global ArrivalNF;
nElenew = nAllEle;
Pnew = Pres;
DD_new = DD;
AllElenew = AllEle;
ConnListnew = ConnList;
global TipStates;
disp('Judge which way the HF will continue to grow');
for i = 1 : nAllEle
    if abs(AllEle(i,10)-3)<1e-6 && Pnew(i)< Mat.Pp*0.1 && ArrivalNF(ConnList(i,1)) > 0.1
        Pnew(i) = 0;%Mat.Pp;
    end
end
for i = 1 : nAllEle
    if TipStates(i) == 999
        %说明i单元是压裂裂缝上与天然裂缝接触的单元
        G = zeros(3,1);
        Num = zeros(3,1);
        count  = 1;
        flag = 0;
        for jj = 1 : ConnList(i,2)
            j = ConnList(i,jj+2);
            if AllEle(j,10) > 1.1 && AllEle(j,10) ~= 2
                %说明j单元已经打开了，就肯定会沿着j单元走
                flag = 1;
                nElenew = nAllEle;
                AllElenew = AllEle;
                ConnListnew = ConnList;
                DD_new = DD;
                Pnew(j) = Pres(i);
                TipStates(i) = 0;
            end
        end
        if flag > 0.1
            return;
        else
            for jj = 1 : ConnList(i,2)
                j = ConnList(i,jj+2);
                if AllEle(j,10) > 1.1
                    Pnew(j) = Pres(i);
                    [DD2,PP2,~]=FluidSolidFullCouple3(EleType,Mat,Pnew(1:nAllEle),DD,dt,nAllEle,AllEle,ConnList,nFracture,Fractures,Sxx,Syy,Sxy,nwell,well,0);
                    Ds= DD2(j);
                    Dn= DD2(j + nAllEle);
                    d = AllEle(j,7)/2;
                    KI1 = -Mat.G*1e3/4/(1-Mat.miu)*(2*pi/d)^0.5*Dn;
                    KI2= -Mat.G*1e3/4/(1-Mat.miu)*(2*pi/d)^0.5*Ds;
                    G(count) = (1-Mat.miu^2)/Mat.E*(KI1^2+KI2^2)/KICnf;
                    if Dn > 1e-8 || PP2(nElenew) < -1e-6
                        %fracture closed
                        G(count) = 0;
                    end
                    Num(count) = j;
                    count = count + 1;
                    Pnew(j) = 0;
                else
                    %开启新的压裂裂缝
                    nElenew = nElenew+ 1;
                    AllElenew(nElenew,1:11) = AllElenew(i,1:11);
                    AllElenew(nElenew,1:2) = AllElenew(i,3:4);
                    AllElenew(nElenew,3:4) = AllElenew(i,3:4)+[AllElenew(i,6),AllElenew(i,5)]*AllElenew(i,7);
                    AllElenew(nElenew,8) = 0.5*(AllElenew(nElenew,1)+AllElenew(nElenew,3));
                    AllElenew(nElenew,9) = 0.5*(AllElenew(nElenew,2)+AllElenew(nElenew,4));
                    ConnListnew = [ConnListnew;ConnList(i,:)];
                    ConnListnew(i,2) = 4;
                    ConnListnew(i,6) = nElenew;
                    ConnListnew(nElenew,2) = 2;
                    ConnListnew(nElenew,3) = i;
                    ConnListnew(nElenew,4) = -1;
                    TipStates(nElenew) = 8;
                    Pnew(nElenew) = Pres(i);
                    DDnew = [DD(1:nAllEle);DD(i)/3;DD(nAllEle+1:2*nAllEle);DD(i+nAllEle)/3];
                    [DD2,PP2,~]=FluidSolidFullCouple3(EleType,Mat,Pnew,DDnew,dt,nElenew,AllElenew,ConnListnew,nFracture,Fractures,Sxx,Syy,Sxy,nwell,well,0);
                    d = AllElenew(nElenew,7)/2;
                    Ds= DD2(nElenew);
                    Dn= DD2(nElenew + nElenew);
                    KI1 = -Mat.G*1e3/4/(1-Mat.miu)*(2*pi/d)^0.5*Dn;
                    KI2= -Mat.G*1e3/4/(1-Mat.miu)*(2*pi/d)^0.5*Ds;
                    G(count) = (1-Mat.miu^2)/Mat.E*(KI1^2+KI2^2)/KIChf;
                    if Dn > 1e-8 || PP2(nElenew) < -1e-6
                        %fracture closed
                        G(count) = 0;
                    end
                    Num(count) = nElenew;
                    count = count + 1;
                end
            end
            [value,index] = max(G);
            if value >1e-6
                if AllEle(Num(index),10) > 1.1
                    %沿着天然裂缝
                    nElenew = nAllEle;
                    AllElenew = AllEle;
                    ConnListnew = ConnList;
                    DD_new = DD;
                    Pnew = Pres;
                    Pnew(Num(index)) = Pres(i);
                end
                TipStates(i) = 0;
            end
        end
        %穿过天然裂缝
    end
end
end
