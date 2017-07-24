function DrawN(nt,curT)
global  AllEle_global sigmaN_global  
global PicScale isMechActive_global PresF_global;
global  nAct  IndexInv    
global FILEPATH;
figure(28);
hold off;
AllEle = AllEle_global;
for i = 1 : nAct
    Dn = sigmaN_global(IndexInv(i))/1e6;%(sigmaN_global(IndexInv(i))-PresF_global(IndexInv(i))*1e6)/1e6;
  %  plot([AllEle(IndexInv(i),1) AllEle(IndexInv(i),3)],[AllEle(IndexInv(i),2) AllEle(IndexInv(i),4)],'c','Linewidth',0.5);
    hold on;
%     if isMechActive_global((i)) == 0
%         plot(AllEle(IndexInv(i),8) ,AllEle(IndexInv(i),9),'m.','Markersize',11);
%     end
%     if isMechActive_global((i)) == -2
%         plot(AllEle(IndexInv(i),8) ,AllEle(IndexInv(i),9),'k.','Markersize',11);
%     end
%     if isMechActive_global((i)) == 1
%         plot(AllEle(IndexInv(i),8) ,AllEle(IndexInv(i),9),'b.','Markersize',11);
%     end
    if abs(sigmaN_global(IndexInv(i))) < 1
        continue;
    end
    scatter(AllEle(IndexInv(i),8),AllEle(IndexInv(i),9),20,Dn,'fill');
%     if Dn > 1e-10
%         plot([AllEle(IndexInv(i),8)+abs(Dn)/2*AllEle(IndexInv(i),5),AllEle(IndexInv(i),8)-abs(Dn)/2*AllEle(IndexInv(i),5)],[AllEle(IndexInv(i),9)-abs(Dn)/2*AllEle(IndexInv(i),6),AllEle(IndexInv(i),9)+abs(Dn)/2*AllEle(IndexInv(i),6)],'b','LineWidth',1);
%     else
%         plot([AllEle(IndexInv(i),8)+abs(Dn)/2*AllEle(IndexInv(i),5),AllEle(IndexInv(i),8)-abs(Dn)/2*AllEle(IndexInv(i),5)],[AllEle(IndexInv(i),9)-abs(Dn)/2*AllEle(IndexInv(i),6),AllEle(IndexInv(i),9)+abs(Dn)/2*AllEle(IndexInv(i),6)],'r','LineWidth',1);
%     end
end

% axis([-0.03 0.06 0.14 0.21]);
axis equal;
%axis(PicScale);
% axis([-0.6 0.6 2.5 4]);
title(['Opening @ Time ',num2str(curT)],'Fontsize',14);
Lx = PicScale(2)-PicScale(1);
Ly = PicScale(4)-PicScale(3);
n = (Ly/Lx);
saveas(28,[FILEPATH,num2str(nt),'sigmaN.fig']);
%close(1);
%close(3);
end