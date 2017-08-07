function DrawE(nt,curT)
global nAllEle_global AllEle_global DD_global nwell well
global PicScale isMechActive_global;
global Index  nAct MaxEle IndexInv  e_global ConnList_global ArrivalNF
global FILEPATH;
figure(20);
hold off;
AllEle = AllEle_global;
nAllEle = nAllEle_global;
Ds = zeros(nAct,1);
Dn = zeros(nAct,1);
for i = 1 : nAct
    Ds(i) = DD_global(IndexInv(i));
    Dn(i) =  e_global(IndexInv(i));%
end
Ds = Ds*2000;
Dn = Dn*8000;
color = ['b','r','g','k'];%Type 1 2 3 4
for i = 1 : nAllEle
    plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b','Linewidth',0.5);
    hold on;
    plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'r.','Markersize',4);
    iparent = ConnList_global(i,1);
    if AllEle(i,10) > 1.1
        if ArrivalNF(iparent) > 0.1
            plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'k','Linewidth',1.5);
            hold on
        end
    end
    if (abs(AllEle(i,10) - 3) < 1e-6) ||  abs(AllEle(i,10) - 1) < 1e-6 || abs(AllEle(i,10) - 6) < 1e-6
        if Index(i) < 0.1
            continue;
        end
         if isMechActive_global(Index(i)) == 0
            plot(AllEle(i,8) ,AllEle(i,9),'m.','Markersize',11);
        end
        if isMechActive_global(Index(i)) == -2
            plot(AllEle(i,8) ,AllEle(i,9),'k.','Markersize',11);
        end
        if isMechActive_global(Index(i)) == 1
            plot(AllEle(i,8) ,AllEle(i,9),'b.','Markersize',11);
        end
        if Dn(Index(i)) < -1e-10
            plot([AllEle(i,8)+abs(Dn(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(Dn(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(Dn(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(Dn(Index(i)))/2*AllEle(i,6)],'b','LineWidth',1);
        else
            plot([AllEle(i,8)+abs(Dn(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(Dn(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(Dn(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(Dn(Index(i)))/2*AllEle(i,6)],'r','LineWidth',1.0);
        end
    end
end
for i = 1 : nwell
         plot([-5   5],[   0 0],'k','Linewidth',2);

    plot([well{i}.heel(1),well{i}.toe(1)],[well{i}.heel(2),well{i}.toe(2)],'k','Linewidth',2);
    for j = 1 : well{i}.nPerf
        num = well{i}.Perfindex(j);
        plot([AllEle(num,1) AllEle(num,3)],[AllEle(num,2) AllEle(num,4)],'r','Linewidth',3);
        plot(well{i}.Perf(j,1),well{i}.Perf(j,2),'wo');
        plot(well{i}.Perf(j,1),well{i}.Perf(j,2),'g.','Markersize',20);
        hold on;
    end
end
% axis([-0.03 0.06 0.14 0.21]);
axis equal;
axis(PicScale);
% axis([-0.6 0.6 2.5 4]);
title(['Opening @ Time ',num2str(curT)],'Fontsize',14);
Lx = PicScale(2)-PicScale(1);
Ly = PicScale(4)-PicScale(3);
n = (Ly/Lx);
set(18,'Position',[20 100 500 500*n])
saveas(18,[FILEPATH,num2str(nt),'E.fig']);
%close(1);
%close(3);
end