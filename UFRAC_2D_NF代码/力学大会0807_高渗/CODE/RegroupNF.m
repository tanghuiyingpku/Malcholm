function  [isRegroup,AllElenew] = RegroupNF(Fractures,nEle,AllEle,ConnList,DD,P)
%Check天然裂缝的状态 是开启、slide还是sticking
% DD is of size 2*nALLELE
global ArrivalNF;
AllElenew = AllEle;
DDs = DD(1:nEle);
DDn = DD(nEle+1: nEle*2);
ElementType = AllEle(1:nEle,10);
for ii = 1 : nEle
    if AllEle(ii,10) > 1.1 %除了压裂裂缝之外的形态
        parent = ConnList(ii,1);
        if ArrivalNF(parent) > 0.1 && AllEle(ii,10) == 2
            continue;
        end
        Dn = DDn(ii);
        Ds = DDs(ii);
        Ks = Fractures{parent}.Ks;
        Kn = Fractures{parent}.Kn;
        coh = Fractures{parent}.Coh;
        frc = Fractures{parent}.Fang;
        if Dn < -1e-16
            % Judge for Dn of NF, if Dn < 0 --> open
            AllElenew(ii,10) = 3;
            AllElenew(ii,11) = 0;
        end
        if Dn > 1e-16
            if abs(AllEle(ii,10) - 2)<1e-6
                %From stick to slip
                ts = -Ks*Ds;
                tn = -Kn*Dn;
                if P(ii) < 1e-6
                    P(ii) = 0;
                end
                if abs(ts) < coh + (abs(tn)-P(ii))*tand(frc)
                    % not slipe stick
                    AllElenew(ii,10) = 2;
                else
                    % slipe
                    AllElenew(ii,10) = 4;
                end
                AllElenew(ii,11) = ts;
            end
            if abs(AllEle(ii,10) - 3)<1e-6
                %From stick to slip
                ts = -Ks*Ds;
                tn = -Kn*Dn;
                if P(ii) < 1e-6
                    P(ii) = 0;
                end
                if abs(ts) < coh + (abs(tn)-P(ii))*tand(frc)
                    % not slipe stick
                    AllElenew(ii,10) = 2;
                else
                    % slipe
                    AllElenew(ii,10) = 4;
                end
                AllElenew(ii,11) = ts;
            end
            if abs(AllEle(ii,10) - 4)<1e-6
                % slipe
                tn = - Kn*Dn;
                ts =  -Ks*Ds;
                if P(ii) < 1e-6
                    P(ii) = 0;
                end
                if abs(ts) < coh + (abs(tn)-P(ii))*tand(frc)
                    % not slipe stick
                    AllElenew(ii,10) = 2;
                end
            end
        end
    else
%         if DDn(ii) > 1e-8
%         %压裂裂缝关闭
%             AllElenew(ii,10) = 2;
%         end
    end
end
EleType_new = AllElenew(1:nEle,10);
if norm(ElementType - EleType_new) > 0.1
    isRegroup = 1;
else
    isRegroup = 0;
end
end
