function DrawDisplacement_new()
%(nEle,DD,AllEle,ConnList,Tipcoord,d,nf)
global  Index IndexInv FILEPATH nAct DD_global nAllEle_global;
global  Mat PicScale AllEle_global ConnList_global MaxEle;

nEle = nAct;
DD = zeros(nEle*2,1);
AllEle = zeros(nEle,11);
ConnList = zeros(nEle,8);

for i = 1 : nAct
    DD(i) = DD_global(IndexInv(i));
    DD(i+nEle) = DD_global(MaxEle+IndexInv(i));
    AllEle(i,:)  = AllEle_global(IndexInv(i),:);
    ConnList(i,:) = ConnList_global(IndexInv(i),:);
end
disp('Show Displacement and Stress around Fracture Surface');
disp('');
eps = 1e-4;

xx = zeros(nEle*2,1);
yy = zeros(nEle*2,1);
figure(44);
hold off;
for i = 1 : nEle
    xx((i-1)*2+1) = AllEle((i),8)-eps*AllEle((i),5);
    xx(i*2) = AllEle((i),8)+eps*AllEle((i),5);
    yy((i-1)*2+1) = AllEle((i),9)+eps*AllEle((i),6);
    yy(i*2) = AllEle((i),9)-eps*AllEle((i),6);
end
[Ux,Uy,~,~,~,~,~] = CalcPointStress_C_local(Mat,DD,xx,yy,AllEle);

amp = 600; %;0.8/max(abs(Ux)+abs(Uy));
% for i = 1 : nEle*2
%     % for j = 1 : nEle*2
%     quiver(xx(i),yy(i),Ux(i)*amp,Uy(i)*amp,'k');
%     hold on;
%     plot(xx(i),yy(i),'k.');
%     % end
% end

for i = 1 : nEle
    nb = ConnList(i,2);
    numi = (i-1)*2+1;
%     if IndexInv(i) > 0.1
%     if TipStates(IndexInv(i)) > 0.1 && TipStates(IndexInv(i)) < 999
%         plot([xx(numi)+Ux(numi)*amp,Tipcoordinate(IndexInv(i),1)],[yy(numi)+Uy(numi)*amp,Tipcoordinate(IndexInv(i),2)],'k','Linewidth',2);
%         hold on;
%         plot([xx(numi)+Ux(numi)*amp,Tipcoordinate(IndexInv(i),1)],[yy(numi)+Uy(numi)*amp,Tipcoordinate(IndexInv(i),2)],'r.');
%         plot([xx(numi+1)+Ux(numi+1)*amp,Tipcoordinate(IndexInv(i),1)],[yy(numi+1)+Uy(numi+1)*amp,Tipcoordinate(IndexInv(i),2)],'k','Linewidth',2);
%         plot([xx(numi+1)+Ux(numi+1)*amp,Tipcoordinate(IndexInv(i),1)],[yy(numi+1)+Uy(numi+1)*amp,Tipcoordinate(IndexInv(i),2)],'r.');
%     end
%     end
    for j = 1 : nb
        if ConnList(i,2+j) < 0.1
            continue;
        end
        if Index(ConnList(i,2+j)) > 0.1
            numj = Index(ConnList(i,2+j));
            if numj > 998
                numj = nEle;
            end
            numj = (numj-1)*2+1;
            if ConnList(i,2+j) < 999
                
                dis1 = (xx(numi) - xx(numj))^2 + (yy(numi) - yy(numj))^2;
                dis2 = (xx(numi) - xx(numj+1))^2 + (yy(numi) - yy(numj+1))^2;
                
                if dis1 < dis2
                   % plot([xx(numi)+Ux(numi)*amp,xx(numj)+Ux(numj)*amp],[yy(numi)+Uy(numi)*amp,yy(numj)+Uy(numj)*amp],'k','Linewidth',2);
                    plot([xx(numi)+Ux(numi)*amp,xx(numj)+Ux(numj)*amp],[yy(numi)+Uy(numi)*amp,yy(numj)+Uy(numj)*amp],'r.','Markersize',11);
                    hold on;
                    %plot([xx(numi+1)+Ux(numi+1)*amp,xx(numj+1)+Ux(numj+1)*amp],[yy(numi+1)+Uy(numi+1)*amp,yy(numj+1)+Uy(numj+1)*amp],'k','Linewidth',2);
                    plot([xx(numi+1)+Ux(numi+1)*amp,xx(numj+1)+Ux(numj+1)*amp],[yy(numi+1)+Uy(numi+1)*amp,yy(numj+1)+Uy(numj+1)*amp],'r.','Markersize',11);
                else
                  %  plot([xx(numi)+Ux(numi)*amp,xx(numj+1)+Ux(numj+1)*amp],[yy(numi)+Uy(numi)*amp,yy(numj+1)+Uy(numj+1)*amp],'k','Linewidth',2);
                    hold on
                    plot([xx(numi)+Ux(numi)*amp,xx(numj+1)+Ux(numj+1)*amp],[yy(numi)+Uy(numi)*amp,yy(numj+1)+Uy(numj+1)*amp],'r.','Markersize',11);
                    
                   % plot([xx(numi+1)+Ux(numi+1)*amp,xx(numj)+Ux(numj)*amp],[yy(numi+1)+Uy(numi+1)*amp,yy(numj)+Uy(numj)*amp],'k','Linewidth',2);
                    plot([xx(numi+1)+Ux(numi+1)*amp,xx(numj)+Ux(numj)*amp],[yy(numi+1)+Uy(numi+1)*amp,yy(numj)+Uy(numj)*amp],'r.','Markersize',11);
                end
                
            else
                numj = (nEle-1)*2+1;
                dis1 = (xx(numi) - xx(numj))^2 + (yy(numi) - yy(numj))^2;
                dis2 = (xx(numi) - xx(numj+1))^2 + (yy(numi) - yy(numj+1))^2;
                
                if dis1 < dis2
                    plot([xx(numi)+Ux(numi)*amp,xx(numj)+Ux(numj)*amp],[yy(numi)+Uy(numi)*amp,yy(numj)+Uy(numj)*amp],'k.','Linewidth',2);
                    hold on
                    plot([xx(numi)+Ux(numi)*amp,xx(numj)+Ux(numj)*amp],[yy(numi)+Uy(numi)*amp,yy(numj)+Uy(numj)*amp],'r.');
                    
                    plot([xx(numi+1)+Ux(numi+1)*amp,xx(numj+1)+Ux(numj+1)*amp],[yy(numi+1)+Uy(numi+1)*amp,yy(numj+1)+Uy(numj+1)*amp],'k.','Linewidth',2);
                    plot([xx(numi+1)+Ux(numi+1)*amp,xx(numj+1)+Ux(numj+1)*amp],[yy(numi+1)+Uy(numi+1)*amp,yy(numj+1)+Uy(numj+1)*amp],'r.');
                else
                    plot([xx(numi)+Ux(numi)*amp,xx(numj+1)+Ux(numj+1)*amp],[yy(numi)+Uy(numi)*amp,yy(numj+1)+Uy(numj+1)*amp],'k.','Linewidth',2);
                    hold on
                    plot([xx(numi)+Ux(numi)*amp,xx(numj+1)+Ux(numj+1)*amp],[yy(numi)+Uy(numi)*amp,yy(numj+1)+Uy(numj+1)*amp],'r.');
                    
                    plot([xx(numi+1)+Ux(numi+1)*amp,xx(numj)+Ux(numj)*amp],[yy(numi+1)+Uy(numi+1)*amp,yy(numj)+Uy(numj)*amp],'k.','Linewidth',2);
                    plot([xx(numi+1)+Ux(numi+1)*amp,xx(numj)+Ux(numj)*amp],[yy(numi+1)+Uy(numi+1)*amp,yy(numj)+Uy(numj)*amp],'r.');
                end
                    
            end
        end
    end
%     if i == nEle
%         epst = d/2;
%         if abs(AllEle(i,1) - Tipcoord(1)) > 1e-6
%             [Uxi,Uyi,~,~,~,~,~] =  CalcPointStress_C_local(Mat,DD,AllEle(i,1)+epst*AllEle(i,6),AllEle(i,2)++epst*AllEle(i,5),AllEle);
%             plot([xx(numi)+Ux(numi)*amp,AllEle(i,1)+epst*AllEle(i,6)+Uxi*amp],[yy(numi)+Uy(numi)*amp,AllEle(i,2)++epst*AllEle(i,5)+Uyi*amp],'k','Linewidth',2);
%             plot([xx(numi)+Ux(numi)*amp,AllEle(i,1)+epst*AllEle(i,6)+Uxi*amp],[yy(numi)+Uy(numi)*amp,AllEle(i,2)++epst*AllEle(i,5)+Uyi*amp],'r.');
%             plot([xx(numi+1)+Ux(numi+1)*amp,AllEle(i,1)+epst*AllEle(i,6)+Uxi*amp],[yy(numi+1)+Uy(numi+1)*amp,AllEle(i,2)+epst*AllEle(i,5)+Uyi*amp],'k','Linewidth',2);
%             plot([xx(numi+1)+Ux(numi+1)*amp,AllEle(i,1)+epst*AllEle(i,6)+Uxi*amp],[yy(numi+1)+Uy(numi+1)*amp,AllEle(i,2)+epst*AllEle(i,5)+Uyi*amp],'r.');
%         else
%             [Uxi,Uyi,~,~,~,~,~] =  CalcPointStress_C_local(Mat,DD,AllEle(i,3)+epst*AllEle(i,6),AllEle(i,4)++epst*AllEle(i,5),AllEle);
%             plot([xx(numi)+Ux(numi)*amp,AllEle(i,3)+epst*AllEle(i,6)+Uxi*amp],[yy(numi)+Uy(numi)*amp,AllEle(i,4)+epst*AllEle(i,5)+Uyi*amp],'k','Linewidth',2);
%             plot([xx(numi)+Ux(numi)*amp,AllEle(i,3)+epst*AllEle(i,6)+Uxi*amp],[yy(numi)+Uy(numi)*amp,AllEle(i,4)+epst*AllEle(i,5)+Uyi*amp],'r.');
%             plot([xx(numi+1)+Ux(numi+1)*amp,AllEle(i,3)+epst*AllEle(i,6)+Uxi*amp],[yy(numi+1)+Uy(numi+1)*amp,AllEle(i,4)+epst*AllEle(i,5)+Uyi*amp],'k','Linewidth',2);
%             plot([xx(numi+1)+Ux(numi+1)*amp,AllEle(i,3)+epst*AllEle(i,6)+Uxi*amp],[yy(numi+1)+Uy(numi+1)*amp,AllEle(i,4)+epst*AllEle(i,5)+Uyi*amp],'r.');
%         end
%     end
end
% Draw NF
num = 0;
xx2 = zeros(500,1);
yy2 = zeros(500,1);
for i = 1 : nAllEle_global
    if AllEle_global(i,10) > 1.1 && Index(i) < 0.1
        num = num+ 1;
        xx2(num) = AllEle_global(i,8);
        yy2(num) = AllEle_global(i,9);
    end
end
xx2 = xx2(1:num);
yy2 = yy2(1:num);
[a,b] = sort(xx2);
xx2 = xx2(b);
yy2 = yy2(b);
[Ux2,Uy2,~,~,~,~,~] =  CalcPointStress_C_local(Mat,DD,xx2,yy2,AllEle);
for numi = 1 : num-1
     plot([xx2(numi)+Ux2(numi)*amp,xx2(numi+1)+Ux2(numi+1)*amp],[yy2(numi)+Uy2(numi)*amp,yy2(numi+1)+Uy2(numi+1)*amp],'k.','Linewidth',2);
end
axis equal;
axis(PicScale);
set(gca,'fontsize',14)
title('Displacement at intersection','Fontsize',13);
saveas(44,[FILEPATH,'Disp600.fig']);

end
