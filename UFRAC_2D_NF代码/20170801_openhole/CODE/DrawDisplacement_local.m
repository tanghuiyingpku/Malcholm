function DrawDisplacement_local(nEle,DD,AllEle,ConnList,Tipcoord,d,nf)
global  TipStates Tipcoordinate  ;
global  Index IndexInv;
global  Mat nAllEle_global AllEle_global ;

disp('Show Displacements of fractures at intersection trials');
disp('');
eps = 1e-12;
xx = zeros(nEle*2,1);
yy = zeros(nEle*2,1);

figure(33);

for i = 1 : nEle
    xx((i-1)*2+1) = AllEle(i,8)-eps*AllEle(i,5);
    xx(i*2) = AllEle(i,8)+eps*AllEle(i,5);
    yy((i-1)*2+1) = AllEle(i,9)+eps*AllEle(i,6);
    yy(i*2) = AllEle(i,9)-eps*AllEle(i,6);
end
[Ux,Uy,~,~,~,~,~] =  CalcPointStress_C_local(Mat,DD,xx,yy,AllEle);

amp = 0.1/max(abs(Ux)+abs(Uy));
%amp = 1;
% Draw NF
num = 0;
xx2 = zeros(500,1);
yy2 = zeros(500,1);
for i = 1 : nAllEle_global
    if AllEle_global(i,10) > 1.1 && i~= nf
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
     plot([xx2(numi)+Ux2(numi)*amp,xx2(numi+1)+Ux2(numi+1)*amp],[yy2(numi)+Uy2(numi)*amp,yy2(numi+1)+Uy2(numi+1)*amp],'m','Linewidth',2);
end

for i = 1 : nEle
    nb = ConnList(i,2);
    numi = (i-1)*2+1;
    if IndexInv(i) > 0.1
    if TipStates(IndexInv(i)) > 0.1 && TipStates(IndexInv(i)) < 999
        plot([xx(numi)+Ux(numi)*amp,Tipcoordinate(IndexInv(i),1)],[yy(numi)+Uy(numi)*amp,Tipcoordinate(IndexInv(i),2)],'k','Linewidth',2);
        hold on;
        plot([xx(numi)+Ux(numi)*amp,Tipcoordinate(IndexInv(i),1)],[yy(numi)+Uy(numi)*amp,Tipcoordinate(IndexInv(i),2)],'r.');
        plot([xx(numi+1)+Ux(numi+1)*amp,Tipcoordinate(IndexInv(i),1)],[yy(numi+1)+Uy(numi+1)*amp,Tipcoordinate(IndexInv(i),2)],'k','Linewidth',2);
        plot([xx(numi+1)+Ux(numi+1)*amp,Tipcoordinate(IndexInv(i),1)],[yy(numi+1)+Uy(numi+1)*amp,Tipcoordinate(IndexInv(i),2)],'r.');
    end
    end
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
                plot([xx(numi)+Ux(numi)*amp,xx(numj)+Ux(numj)*amp],[yy(numi)+Uy(numi)*amp,yy(numj)+Uy(numj)*amp],'k','Linewidth',2);
                hold on;
                plot([xx(numi)+Ux(numi)*amp,xx(numj)+Ux(numj)*amp],[yy(numi)+Uy(numi)*amp,yy(numj)+Uy(numj)*amp],'r.');
                plot([xx(numi+1)+Ux(numi+1)*amp,xx(numj+1)+Ux(numj+1)*amp],[yy(numi+1)+Uy(numi+1)*amp,yy(numj+1)+Uy(numj+1)*amp],'k','Linewidth',2);
                plot([xx(numi+1)+Ux(numi+1)*amp,xx(numj+1)+Ux(numj+1)*amp],[yy(numi+1)+Uy(numi+1)*amp,yy(numj+1)+Uy(numj+1)*amp],'r.');
            else
                numj = (nEle-1)*2+1;
                plot([xx(numi)+Ux(numi)*amp,xx(numj)+Ux(numj)*amp],[yy(numi)+Uy(numi)*amp,yy(numj)+Uy(numj)*amp],'k','Linewidth',2);
                hold on
                plot([xx(numi)+Ux(numi)*amp,xx(numj)+Ux(numj)*amp],[yy(numi)+Uy(numi)*amp,yy(numj)+Uy(numj)*amp],'r.');
                plot([xx(numi+1)+Ux(numi+1)*amp,xx(numj+1)+Ux(numj+1)*amp],[yy(numi+1)+Uy(numi+1)*amp,yy(numj+1)+Uy(numj+1)*amp],'k','Linewidth',2);
                plot([xx(numi+1)+Ux(numi+1)*amp,xx(numj+1)+Ux(numj+1)*amp],[yy(numi+1)+Uy(numi+1)*amp,yy(numj+1)+Uy(numj+1)*amp],'r.');
            end
        end
    end
    if i == nEle
        epst = d/2;
        if abs(AllEle(i,1) - Tipcoord(1)) > 1e-6
            [Uxi,Uyi,~,~,~,~,~] =  CalcPointStress_C_local(Mat,DD,AllEle(i,1)+epst*AllEle(i,6),AllEle(i,2)++epst*AllEle(i,5),AllEle);
            plot([xx(numi)+Ux(numi)*amp,AllEle(i,1)+epst*AllEle(i,6)+Uxi*amp],[yy(numi)+Uy(numi)*amp,AllEle(i,2)++epst*AllEle(i,5)+Uyi*amp],'k','Linewidth',2);
            plot([xx(numi)+Ux(numi)*amp,AllEle(i,1)+epst*AllEle(i,6)+Uxi*amp],[yy(numi)+Uy(numi)*amp,AllEle(i,2)++epst*AllEle(i,5)+Uyi*amp],'r.');
            plot([xx(numi+1)+Ux(numi+1)*amp,AllEle(i,1)+epst*AllEle(i,6)+Uxi*amp],[yy(numi+1)+Uy(numi+1)*amp,AllEle(i,2)+epst*AllEle(i,5)+Uyi*amp],'k','Linewidth',2);
            plot([xx(numi+1)+Ux(numi+1)*amp,AllEle(i,1)+epst*AllEle(i,6)+Uxi*amp],[yy(numi+1)+Uy(numi+1)*amp,AllEle(i,2)+epst*AllEle(i,5)+Uyi*amp],'r.');
        else
            [Uxi,Uyi,~,~,~,~,~] =  CalcPointStress_C_local(Mat,DD,AllEle(i,3)+epst*AllEle(i,6),AllEle(i,4)++epst*AllEle(i,5),AllEle);
            plot([xx(numi)+Ux(numi)*amp,AllEle(i,3)+epst*AllEle(i,6)+Uxi*amp],[yy(numi)+Uy(numi)*amp,AllEle(i,4)+epst*AllEle(i,5)+Uyi*amp],'k','Linewidth',2);
            plot([xx(numi)+Ux(numi)*amp,AllEle(i,3)+epst*AllEle(i,6)+Uxi*amp],[yy(numi)+Uy(numi)*amp,AllEle(i,4)+epst*AllEle(i,5)+Uyi*amp],'r.');
            plot([xx(numi+1)+Ux(numi+1)*amp,AllEle(i,3)+epst*AllEle(i,6)+Uxi*amp],[yy(numi+1)+Uy(numi+1)*amp,AllEle(i,4)+epst*AllEle(i,5)+Uyi*amp],'k','Linewidth',2);
            plot([xx(numi+1)+Ux(numi+1)*amp,AllEle(i,3)+epst*AllEle(i,6)+Uxi*amp],[yy(numi+1)+Uy(numi+1)*amp,AllEle(i,4)+epst*AllEle(i,5)+Uyi*amp],'r.');
        end
    end
end

axis equal;
axis([-0.03 0.03 0.45 0.49]);
title('Displacement at intersection','Fontsize',13);

end
