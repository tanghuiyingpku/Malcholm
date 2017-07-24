function DrawPressure(nt,curT)
global nAllEle_global AllEle_global  nwell well  ConnList_global ArrivalNF
global PicScale;
global Index  nAct  IndexInv PresF_global CpF_global
global FILEPATH Mat;

AllEle = AllEle_global;
nAllEle = nAllEle_global;
PP = zeros(nAct,1);
Cp = zeros(nAct,1);
for i = 1 : nAct
    PP(i) = PresF_global(IndexInv(i));
    Cp(i) = CpF_global(IndexInv(i),1);
end
figure(2);
hold off
%PP = (PP)+Mat.Sxx;
PP =PP/2;
for i = 1 : nAllEle
     plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b','Linewidth',1.5);
    hold on;
    plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'r.','Linewidth',1.5);
    iparent = ConnList_global(i,1);
    if AllEle(i,10) > 1.1
        if ArrivalNF(iparent) > 0.1
            plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'k','Linewidth',1.5);
            hold on
        end
    end
    if (abs(AllEle(i,10) - 3) < 1e-6) ||  abs(AllEle(i,10) - 1) < 1e-6
        if Index(i) < 0.1
            continue;
        end
        plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b','Linewidth',1.5);
        hold on;
        if PP(Index(i)) < -1e-3
            plot([AllEle(i,8)+abs(PP(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(PP(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(PP(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(PP(Index(i)))/2*AllEle(i,6)],'k','LineWidth',1.5);
        else
            plot([AllEle(i,8)+abs(PP(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(PP(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(PP(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(PP(Index(i)))/2*AllEle(i,6)],'m','LineWidth',1.5);
        end
        hold on;
       % plot([AllEle(i,8)+abs(PP(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(PP(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(PP(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(PP(Index(i)))/2*AllEle(i,6)],'b.','Markersize',12);
    end
    hold on;
%     plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b');
%     plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'k.');
%     plot(AllEle(i,8),AllEle(i,9),'r*');
end
for i = 1 : nwell
    plot([well{i}.heel(1),well{i}.toe(1)],[well{i}.heel(2),well{i}.toe(2)],'k','Linewidth',2);
    for j = 1 : well{i}.nPerf
        num = well{i}.Perfindex(j);
        plot([AllEle(num,1) AllEle(num,3)],[AllEle(num,2) AllEle(num,4)],'r','Linewidth',3);
        plot(well{i}.Perf(j,1),well{i}.Perf(j,2),'wo');
        plot(well{i}.Perf(j,1),well{i}.Perf(j,2),'g.','Markersize',20);
        hold on;
    end
end
title(['Presure@ Time ',num2str(curT)],'Fontsize',14);
axis equal;
axis(PicScale);
Lx = PicScale(2)-PicScale(1);
Ly = PicScale(4)-PicScale(3);
n = (Ly/Lx);
saveas(2,[FILEPATH,num2str(nt),'Pres.fig']);
%close(2);
% figure(6);
% hold off
% Cp = Cp * 50;
% for i = 1 : nAllEle
%     if (abs(AllEle(i,10) - 3) < 1e-6) ||  abs(AllEle(i,10) - 1) < 1e-6  && Index(i) > 0.1
%         if Index(i) < 0.1
%             f = 1;
%         end
%         if Cp(Index(i)) < -1e-3
%             plot([AllEle(i,8)+abs(Cp(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(Cp(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(Cp(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(Cp(Index(i)))/2*AllEle(i,6)],'k');
%         else
%             plot([AllEle(i,8)+abs(Cp(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(Cp(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(Cp(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(Cp(Index(i)))/2*AllEle(i,6)],'m');
%         end
%         hold on;
%         plot([AllEle(i,8)+abs(Cp(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(Cp(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(Cp(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(Cp(Index(i)))/2*AllEle(i,6)],'b.');
%     end
%     hold on;
%     plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b');
%     plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'k.');
%     plot(AllEle(i,8),AllEle(i,9),'r*');
% end
% for i = 1 : nwell
%     plot([well{i}.heel(1),well{i}.toe(1)],[well{i}.heel(2),well{i}.toe(2)],'k','Linewidth',2);
%     for j = 1 : well{i}.nPerf
%         num = well{i}.Perfindex(j);
%         plot([AllEle(num,1) AllEle(num,3)],[AllEle(num,2) AllEle(num,4)],'r','Linewidth',3);
%         plot(well{i}.Perf(j,1),well{i}.Perf(j,2),'wo');
%         plot(well{i}.Perf(j,1),well{i}.Perf(j,2),'g.','Markersize',20);
%         hold on;
%         
%     end
% end
% title(['Proppant Concentration @ Time ',num2str(curT)],'Fontsize',14);
% axis equal;
% axis(PicScale);
% Lx = PicScale(2)-PicScale(1);
% Ly = PicScale(4)-PicScale(3);
% n = (Ly/Lx);
% set(6,'Position',[300 300 600 600*n])
% saveas(6,[FILEPATH,num2str(nt),'CP.fig']);
%close(6);
end