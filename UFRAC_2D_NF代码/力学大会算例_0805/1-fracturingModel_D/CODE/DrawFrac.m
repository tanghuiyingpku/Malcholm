function DrawFrac(nAllEle,AllEle,nwell,well)
global PicScale;

figure(1);
hold off;
for i = 1 : nAllEle
    plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b');
    hold on;
    plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b*');
    plot(AllEle(i,8) , AllEle(i,9),'r.');
    axis equal;
end
for i = 1 : nwell
    plot(well{i}.Loc(1),well{i}.Loc(2),'ro');
    hold on;
    num = well{i}.index;
    plot([AllEle(num,1) AllEle(num,3)],[AllEle(num,2) AllEle(num,4)],'r','Linewidth',3);
end
axis(PicScale);
end