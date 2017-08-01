function DrawPressure_t(nt,curT)
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
PP = (PP);%+Mat.Sxx;
% PP =PP/10;
for i = 1 : nAllEle
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
        plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b','Linewidth',0.8);
        hold on;
        if PP(Index(i)) < -1e-3
            X = AllEle(i,8);
            Y = AllEle(i,9);
            Z =PP(Index(i));
            scatter(X,Y,60,Z);
            scatter(X,Y,40,Z,'fill');
%             plot([AllEle(i,8)+abs(PP(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(PP(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(PP(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(PP(Index(i)))/2*AllEle(i,6)],'k');
        else
            X = AllEle(i,8);
            Y = AllEle(i,9);
            Z =PP(Index(i));
            scatter(X,Y,60,Z,'fill');
%             plot([AllEle(i,8)+abs(PP(Index(i)))/2*AllEle(i,5),AllEle(i,8)-abs(PP(Index(i)))/2*AllEle(i,5)],[AllEle(i,9)-abs(PP(Index(i)))/2*AllEle(i,6),AllEle(i,9)+abs(PP(Index(i)))/2*AllEle(i,6)],'m');
        end
    end
    hold on;
%     plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b');
%     plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'k.');
%     plot(AllEle(i,8),AllEle(i,9),'r*');
end

title(['Presure@ Time ',num2str(curT)/60,'min'],'Fontsize',14);
axis equal;
% axis(PicScale);
colorbar;
Lx = PicScale(2)-PicScale(1);
Ly = PicScale(4)-PicScale(3);
n = (Ly/Lx);

end