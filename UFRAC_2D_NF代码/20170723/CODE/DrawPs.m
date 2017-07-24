function DrawPs(nt,curT)
global  AllEle_global PresF_global  
global PicScale isMechActive_global;
global   nAct  IndexInv    
global FILEPATH;
figure(18);
hold off;
AllEle = AllEle_global;
Ds = zeros(nAct,1);
Dn = zeros(nAct,1);
for i = 1 : nAct
    Ds(i) = PresF_global(IndexInv(i));
    Ds = Ds/10;%12000;
    Dn = Ds;
    plot([AllEle(IndexInv(i),1) AllEle(IndexInv(i),3)],[AllEle(IndexInv(i),2) AllEle(IndexInv(i),4)],'b','Linewidth',0.5);
    hold on;
    if isMechActive_global((i)) == 0
        plot(AllEle(IndexInv(i),8) ,AllEle(IndexInv(i),9),'m.','Markersize',11);
    end
    if isMechActive_global((i)) == -2
        plot(AllEle(IndexInv(i),8) ,AllEle(IndexInv(i),9),'k.','Markersize',11);
    end
    if isMechActive_global((i)) == 1
        plot(AllEle(IndexInv(i),8) ,AllEle(IndexInv(i),9),'b.','Markersize',11);
    end
    if Dn(i) < -1e-10
        plot([AllEle(IndexInv(i),8)+abs(Dn((i)))/2*AllEle(IndexInv(i),5),AllEle(IndexInv(i),8)-abs(Dn((i)))/2*AllEle(IndexInv(i),5)],[AllEle(IndexInv(i),9)-abs(Dn((i)))/2*AllEle(IndexInv(i),6),AllEle(IndexInv(i),9)+abs(Dn((i)))/2*AllEle(IndexInv(i),6)],'b','LineWidth',1);
    else
        plot([AllEle(IndexInv(i),8)+abs(Dn((i)))/2*AllEle(IndexInv(i),5),AllEle(IndexInv(i),8)-abs(Dn((i)))/2*AllEle(IndexInv(i),5)],[AllEle(IndexInv(i),9)-abs(Dn((i)))/2*AllEle(IndexInv(i),6),AllEle(IndexInv(i),9)+abs(Dn((i)))/2*AllEle(IndexInv(i),6)],'r','LineWidth',1);
    end
end

% axis([-0.03 0.06 0.14 0.21]);
axis equal;
%axis(PicScale);
% axis([-0.6 0.6 2.5 4]);
title(['Opening @ Time ',num2str(curT)],'Fontsize',14);
Lx = PicScale(2)-PicScale(1);
Ly = PicScale(4)-PicScale(3);
n = (Ly/Lx);
set(18,'Position',[20 100 500 500*n])
saveas(18,[FILEPATH,num2str(nt),'Dss.fig']);
%close(1);
%close(3);
end