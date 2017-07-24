function DrawE_t(nt,curT)
global  AllEle_global    
global PicScale;
global   nAct  IndexInv  CpF_global  e_global
global FILEPATH ;

PP = zeros(nAct,1);
Cp = zeros(nAct,1);
for i = 1 : nAct
    PP(i) = e_global(IndexInv(i));
    Cp(i) = CpF_global(IndexInv(i),1);
end
figure(32);
hold off
PP = (PP);%+Mat.Sxx;
% PP =PP/10;
for i = 1 : nAct
    if 1%(abs(AllEle_global(IndexInv(i),10) - 3) < 1e-6) ||  abs(AllEle_global(IndexInv(i),10) - 1) < 1e-6

        plot([AllEle_global(IndexInv(i),1) AllEle_global(IndexInv(i),3)],[AllEle_global(IndexInv(i),2) AllEle_global(IndexInv(i),4)],'b','Linewidth',0.5);
        hold on;
        if PP((i)) < -1e-3
            X = AllEle_global(IndexInv(i),8);
            Y =AllEle_global(IndexInv(i),9);
            Z =PP((i));
            scatter(X,Y,60,Z);
            scatter(X,Y,40,Z,'fill');
%             plot([AllEle(i,8)+abs(PP(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(PP(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(PP(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(PP(Index(i)))/2*AllEle(i,6)],'k');
        else
            X = AllEle_global(IndexInv(i),8);
            Y = AllEle_global(IndexInv(i),9);
            Z =PP((i));
            scatter(X,Y,60,Z,'fill');
%             plot([AllEle(i,8)+abs(PP(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(PP(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(PP(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(PP(Index(i)))/2*AllEle(i,6)],'m');
        end
    end
    hold on;
%     plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b');
%     plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'k.');
%     plot(AllEle(i,8),AllEle(i,9),'r*');
end

title(['Aperture @ Time ',num2str(curT)/60,'min'],'Fontsize',14);
axis equal;
% axis(PicScale);
colorbar;
Lx = PicScale(2)-PicScale(1);
Ly = PicScale(4)-PicScale(3);
n = (Ly/Lx);
set(2,'Position',[300 0 600 600*n])
saveas(2,[FILEPATH,num2str(nt),'Pres.fig']);

end