function DrawDs(nt,curT)
global nAllEle_global AllEle_global DD_global nwell well
global PicScale
global Index  nAct MaxEle IndexInv  ConnList_global ArrivalNF
figure(17);
hold off;
AllEle = AllEle_global;
nAllEle = nAllEle_global;
Ds = zeros(nAct,1);
Dn = zeros(nAct,1);
for i = 1 : nAct
    Ds(i) = DD_global(IndexInv(i));
    Dn(i) = DD_global(IndexInv(i)+MaxEle);% e_global(IndexInv(i));%
end
Ds = Dn*800;
Dn = Dn*5000;
color = ['b','r','g','k'];%Type 1 2 3 4
for i = 1 : nAllEle
    plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b','Linewidth',0.5);
    hold on;
    plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b.','Linewidth',0.5);
%     plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'r.','Markersize',4);
    iparent = ConnList_global(i,1);
%     text(AllEle((i),8) ,AllEle((i),9),num2str((i)));
    if AllEle(i,10) > 1.1
        if ArrivalNF(iparent) > 0.1
            plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'k','Linewidth',0.5);
            hold on
        end
    end
    if 1%(abs(AllEle(i,10) - 3) < 1e-6) ||  abs(AllEle(i,10) - 1) < 1e-6 || abs(AllEle(i,10) - 6) < 1e-6
        if Index(i) < 0.1
            continue;
        end
%         if isMechActive_global(Index(i)) == 0
%             plot(AllEle(i,8) ,AllEle(i,9),'m.','Markersize',11);
%         end
%         if isMechActive_global((i)) == -2
%         if abs(Ds(i)) > 1e-12
%             plot(AllEle(IndexInv(i),8) ,AllEle(IndexInv(i),9),'m.','Markersize',3);
%         else
%             plot(AllEle(IndexInv(i),8) ,AllEle(IndexInv(i),9),'k.','Markersize',3);
%         end
%     end
%         if isMechActive_global(Index(i)) == 1
%             plot(AllEle(i,8) ,AllEle(i,9),'b.','Markersize',11);
%         end
        
        if Ds(Index(i)) < -1e-10
          %  plot([AllEle(i,8)+abs(Ds(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(Ds(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(Ds(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(Ds(Index(i)))/2*AllEle(i,6)],'b','LineWidth',0.5);
            
            plot([AllEle(i,8)+abs(Ds(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(Ds(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(Ds(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(Ds(Index(i)))/2*AllEle(i,6)],'b','Markersize',4,'LineWidth',1.0);
        else
           % plot([AllEle(i,8)+abs(Ds(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(Ds(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(Ds(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(Ds(Index(i)))/2*AllEle(i,6)],'r','LineWidth',0.5);
            
            plot([AllEle(i,8)+abs(Ds(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(Ds(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(Ds(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(Ds(Index(i)))/2*AllEle(i,6)],'r','Markersize',4,'LineWidth',1.0);
        end
        % plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'k^');
    end
end

% axis([-0.03 0.06 0.14 0.21]);
axis equal;
% axis(PicScale);
% axis([-0.6 0.6 2.5 4]);
title(['Opening @ Time ',num2str(curT)],'Fontsize',14);
Lx = PicScale(2)-PicScale(1);
Ly = PicScale(4)-PicScale(3);
n = (Ly/Lx);
% saveas(17,[FILEPATH,num2str(nt),'Ds.fig']);
%close(1);
%close(3);
end