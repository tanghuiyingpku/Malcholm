function DrawPotentialTip(AllEle,nwell,well)
global PicScale;
global TipStates;
global IndexInv nAct

figure(6);
hold off;
for i = 1 : nAct
    plot([AllEle(IndexInv(i),1) AllEle(IndexInv(i),3)],[AllEle(IndexInv(i),2) AllEle(IndexInv(i),4)],'b','Linewidth',1.5);
    hold on;
    %  plot([AllEle(i,8)+abs(DD(i+nAllEle))/2*AllEle(i,5),AllEle(i,8)-abs(DD(i+nAllEle))/2*AllEle(i,5)],[AllEle(i,9)-abs(DD(i+nAllEle))/2*AllEle(i,6),AllEle(i,9)+abs(DD(i+nAllEle))/2*AllEle(i,6)],'b.');
    if abs(TipStates(i)+2) < 1e-3
        plot([AllEle(IndexInv(i),1) AllEle(IndexInv(i),3)],[AllEle(IndexInv(i),2) AllEle(IndexInv(i),4)],'r^');
    end
end
for i = 1 : nwell
    plot(well{i}.Loc(1),well{i}.Loc(2),'ro','Markersize',6);
    hold on;
    num = well{i}.index;
    plot([AllEle(num,1) AllEle(num,3)],[AllEle(num,2) AllEle(num,4)],'r','Linewidth',1.5);
end
axis(PicScale);
title('Potential Tips');
end