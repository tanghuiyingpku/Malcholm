function DrawDn_ta(nt,curT)
global  AllEle_global    nAllEle_global
global PicScale;
global      MaxEle  DD_global
global FILEPATH ;

figure(33);
hold off
% PP =PP/10;
for i = 1 : nAllEle_global
            plot([AllEle_global((i),1) AllEle_global((i),3)],[AllEle_global((i),2) AllEle_global((i),4)],'b','Linewidth',0.5);

    if (abs(AllEle_global((i),10) - 3) < 1e-6) ||  abs(AllEle_global((i),10) - 1) < 1e-6
        
        plot([AllEle_global((i),1) AllEle_global((i),3)],[AllEle_global((i),2) AllEle_global((i),4)],'b','Linewidth',0.5);
        hold on;
        Di=  -DD_global(i+MaxEle) ;
        X = AllEle_global((i),8);
        Y =AllEle_global((i),9);
        Z = Di;
        scatter(X,Y,20,Z,'fill');
        %             plot([AllEle(i,8)+abs(PP(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(PP(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(PP(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(PP(Index(i)))/2*AllEle(i,6)],'k');
        
    end
    hold on;
    %     plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b');
    %     plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'k.');
    %     plot(AllEle(i,8),AllEle(i,9),'r*');
end

title(['Dn_total@ Time ',num2str(curT)/60,'min'],'Fontsize',14);
axis equal;
% axis(PicScale);
colorbar;
Lx = PicScale(2)-PicScale(1);
Ly = PicScale(4)-PicScale(3);
n = (Ly/Lx);
set(33,'Position',[300 0 600 600*n])
saveas(33,[FILEPATH,num2str(nt),'Dn_t.fig']);

end