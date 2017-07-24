function DrawDisplacement_c()
global  AllEle_global ConnList_global
global  PicScale;
global  nAct  IndexInv   MaxEle   DD_global Index
global  Mat;

nEle = nAct;
AllEle = AllEle_global;
ConnList = ConnList_global;
Ds = zeros(nAct,1);
Dn = zeros(nAct,1);
convexhull = zeros(nEle*2,2);
for i = 1 : nAct
    Ds(i) = DD_global(IndexInv(i));
    Dn(i) =DD_global(IndexInv(i)+MaxEle); % e_global(IndexInv(i));%
end
DD = [Ds;Dn];

disp('Show Displacement and Stress around Fracture Surface');
disp('');
eps = 1e-8;
xx = zeros(nEle*2,1);
yy = zeros(nEle*2,1);
i = 0;
for j = 1 : nEle
    if ConnList(IndexInv(j),1) > 1.1
        i = i+1;
        xx(i) = AllEle(IndexInv(j),8)+eps*AllEle(IndexInv(j),5);
        yy(i) = AllEle(IndexInv(j),9)-eps*AllEle(IndexInv(j),6);
    end
end
xx = xx(1:i);
yy = yy(1:i);
xx = zeros(nEle*2,1);
yy = zeros(nEle*2,1);
figure(4);
hold off;
for i = 1 : nEle
    xx((i-1)*2+1) = AllEle(IndexInv(i),8)-eps*AllEle(IndexInv(i),5);
    xx(i*2) = AllEle(IndexInv(i),8)+eps*AllEle(IndexInv(i),5);
    yy((i-1)*2+1) = AllEle(IndexInv(i),9)+eps*AllEle(IndexInv(i),6);
    yy(i*2) = AllEle(IndexInv(i),9)-eps*AllEle(IndexInv(i),6);
end
[Ux,Uy,Sxx,Syy,Sxy,~,~] = CalcPointStress_CD(Mat,DD,xx,yy,AllEle,ConnList);


axis(PicScale);

amp = 0.1/max(abs(Ux)+abs(Uy));
%amp = 1;
for i = 1 : nEle
    plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b','Linewidth',1.5);
    hold on;
end

for i = 1 : nEle*2
    % for j = 1 : nEle*2
    quiver(xx(i),yy(i),Ux(i)*amp,Uy(i)*amp,'k');
    hold on;
    % end
end
axis(PicScale);
title('Displacement','Fontsize',13);

end
