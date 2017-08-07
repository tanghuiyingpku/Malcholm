function DrawOpen_backup(nt,curT)
global nAllEle_global AllEle_global DD_global nwell well e_global 
global PicScale ConnList_global;
global Index  nAct MaxEle IndexInv
global FILEPATH ArrivalNF;
figure(11);
hold off;
AllEle = AllEle_global;
nAllEle = nAllEle_global;
Ds = zeros(nAct,1);
Dn = zeros(nAct,1);
for i = 1 : nAct
    Ds(i) = DD_global(IndexInv(i));
    Dn(i) = DD_global(IndexInv(i)+MaxEle);%DD_global(IndexInv(i)+MaxEle);
end
Ds = Ds*38800;
Dn = Dn*38800;
color = ['b','r','g','k'];%Type 1 2 3 4
for i = 1 : nAllEle
   % plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b','Linewidth',1.5);
   % hold on;
   iparent = ConnList_global(i,1);
   if iparent < 0.1
       f = 2;
   end
   if AllEle(i,10) > 1.1
       if ArrivalNF(iparent) > 0.1
           plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b','Linewidth',0.8);
           hold on
       end
   end
    if (abs(AllEle(i,10) - 3) < 1e-6) ||  abs(AllEle(i,10) - 1) < 1e-6 
        if Index(i) < 0.1
           continue;
        end
        plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b','Linewidth',0.8);
        hold on;
        if Dn(Index(i)) < -1e-10
            plot([AllEle(i,8)+abs(Dn(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(Dn(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(Dn(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(Dn(Index(i)))/2*AllEle(i,6)],'k','LineWidth',0.5);
        else
            plot([AllEle(i,8)+abs(Dn(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(Dn(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(Dn(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(Dn(Index(i)))/2*AllEle(i,6)],'r','LineWidth',0.5);
        end
        hold on;
        if Ds(Index(i)) < -1e-10
            plot([AllEle(i,8)+abs(Ds(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(Ds(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(Ds(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(Ds(Index(i)))/2*AllEle(i,6)],'k.','Markersize',7);
        else
            plot([AllEle(i,8)+abs(Ds(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(Ds(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(Ds(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(Ds(Index(i)))/2*AllEle(i,6)],'r.','Markersize',7);
        end
        hold on;
       % plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'k^');
    end    
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
axis equal;
%axis(PicScale);
axis([200-3 200+3 200-3 200+3]);
title(['Opening @ Time ',num2str(curT/60),'min'],'Fontsize',14);
Lx = PicScale(2)-PicScale(1);
Ly = PicScale(4)-PicScale(3);
n = (Ly/Lx);
set(11,'Position',[20 0 500 500*n])
saveas(11,[FILEPATH,num2str(nt),'Open.fig']);
%close(1);
figure(3);
hold off;
for i = 1 : nAllEle 
  %  plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],color(AllEle(i,10)),'Marker','o','Markersize',5);
  %  hold on;
   % plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],color(AllEle(i,10)),'LineWidth',3);
   iparent = ConnList_global(i,1);
   if AllEle(i,10) > 1.1
       if ArrivalNF(iparent) > 0.1
           plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],color(AllEle(i,10)),'Marker','o','Markersize',5);
           hold on;
           plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],color(AllEle(i,10)),'LineWidth',3);
       end
   end
    if Index(i) > 0.1
        plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b','Linewidth',1.5);
        hold on;
        plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],color(AllEle(i,10)),'Marker','o','Markersize',5);
        if abs(Dn(Index(i))) < 1e-15 ||  Dn(Index(i)) > 1e-15
            plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'m*','LineWidth',2,'Markersize',5);
            hold on;
            if AllEle(i,11) > 0.1
                plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'k--','LineWidth',4);
            else
                plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'m--','LineWidth',4);
            end
        end
    end
    %  plot([AllEle(i,8)+abs(DD(i+nAllEle))/2*AllEle(i,5),AllEle(i,8)-abs(DD(i+nAllEle))/2*AllEle(i,5)],[AllEle(i,9)-abs(DD(i+nAllEle))/2*AllEle(i,6),AllEle(i,9)+abs(DD(i+nAllEle))/2*AllEle(i,6)],'b.');
end
for i = 1 : nwell
    plot([well{i}.heel(1),well{i}.toe(1)],[well{i}.heel(2),well{i}.toe(2)],'k','Linewidth',2);
    for j = 1 : well{i}.nPerf
        num = well{i}.Perfindex(j);
        plot([AllEle(num,1) AllEle(num,3)],[AllEle(num,2) AllEle(num,4)],'r','Linewidth',3);
        plot(well{i}.Perf(j,1),well{i}.Perf(j,2),'wo');
        plot(well{i}.Perf(j,1),well{i}.Perf(j,2),'g.','Markersize',5);
        hold on;
    end
end
% for i = 1 : nwell
%     plot(well{i}.Loc(1),well{i}.Loc(2),'m.','Markersize',29);
%     hold on;
%     num = well{i}.index;
%     plot([AllEle(num,1) AllEle(num,3)],[AllEle(num,2) AllEle(num,4)],'w+','Linewidth',1.5);
% end
axis equal;
axis(PicScale);
title(['Fracture Path & Type Distribution @ Time ',num2str(curT/60),'min'],'Fontsize',14);
Lx = PicScale(2)-PicScale(1);
Ly = PicScale(4)-PicScale(3);
n = (Ly/Lx);
set(3,'Position',[400 500 500 500*n])
saveas(3,[FILEPATH,num2str(nt),'Type.fig']);
%close(3);
end