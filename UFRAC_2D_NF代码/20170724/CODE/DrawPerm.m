function DrawPerm(nt,curT)
global  AllEle_global sigmaN_global
global PicScale isMechActive_global ;
global  nAct  IndexInv  DD_global   e_global
global FILEPATH;
figure(48);
hold off;
AllEle = AllEle_global;
for i = 1 : nAct
    Dn = e_global(IndexInv(i)) + tand(30)*abs(DD_global(IndexInv(i)));
    Dn = Dn^2/12*1e15;
    plot([AllEle(IndexInv(i),1) AllEle(IndexInv(i),3)],[AllEle(IndexInv(i),2) AllEle(IndexInv(i),4)],'c','Linewidth',0.5);
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
    
    X = AllEle_global(IndexInv(i),8);
    Y = AllEle_global(IndexInv(i),9);
    Z =Dn;
    scatter(X,Y,20,Z,'fill');
    %             plot([AllEle(i,8)+abs(PP(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(PP(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(PP(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(PP(Index(i)))/2*AllEle(i,6)],'m');
    
end

% axis([-0.03 0.06 0.14 0.21]);
axis equal;
%axis(PicScale);
% axis([-0.6 0.6 2.5 4]);
title(['Opening @ Time ',num2str(curT)],'Fontsize',14);
Lx = PicScale(2)-PicScale(1);
Ly = PicScale(4)-PicScale(3);
n = (Ly/Lx);
set(48,'Position',[20 100 500 500*n])
saveas(48,[FILEPATH,num2str(nt),'sigmaN.fig']);
%close(1);
%close(3);
end