function DrawDisplacement(nt,curT)
global nAllEle_global AllEle_global  nwell well  ConnList_global 
global PicScale;
global   nAct  IndexInv   MaxEle nFracture Fractures DD_global
global FILEPATH Mat;
global WatchLine;

nEle = nAct;
AllEle = AllEle_global(1:nAct,:);
ConnList = ConnList_global(1:nAct,:);
SigmaX = Mat.Sxx;
SigmaY = Mat.Syy;
SigmaXY = Mat.Sxy;
Ds = zeros(nAct,1);
Dn = zeros(nAct,1);
for i = 1 : nAct
    AllEle(i,:) = AllEle_global(IndexInv(i),:);
    ConnList(i,:) = ConnList_global(IndexInv(i),:);
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
    if ConnList(j,1) > 1.1
        i = i+1;
        xx(i) = AllEle(j,8)+eps*AllEle(j,5);
        yy(i) = AllEle(j,9)-eps*AllEle(j,6);
    end
end
xx = xx(1:i);
yy = yy(1:i);
xx = zeros(nEle*2,1);
yy = zeros(nEle*2,1);
figure(4);
hold off;
for i = 1 : nEle
    xx((i-1)*2+1) = AllEle(i,8)-eps*AllEle(i,5);
    xx(i*2) = AllEle(i,8)+eps*AllEle(i,5);
    yy((i-1)*2+1) = AllEle(i,9)+eps*AllEle(i,6);
    yy(i*2) = AllEle(i,9)-eps*AllEle(i,6);
end
[Ux,Uy,Sxx,Syy,Sxy,~,~] = CalcPointStress_H(Mat,DD,xx,yy,nFracture,Fractures,nEle,AllEle,ConnList);

Sxx = Sxx + SigmaX;
Syy = Syy + SigmaY;
Sxy = Sxy + SigmaXY;
Un = zeros(nEle*2,1);
Us = zeros(nEle*2,1);

Sn = zeros(nEle*2,1);
St = zeros(nEle*2,1);
Tt = zeros(nEle*2,1);
%Stµƒº∆À„…‘Œ¢∏¥‘”µ„
for i = 1: nEle
    cosb = AllEle(i,6);
    sinb = AllEle(i,5);
    Us((i-1)*2+1:i*2) = Ux((i-1)*2+1:i*2) *cosb + Uy((i-1)*2+1:i*2)*sinb;
    Un((i-1)*2+1:i*2) = -Ux((i-1)*2+1:i*2) *sinb + Uy((i-1)*2+1:i*2)*cosb;
    Sn((i-1)*2+1) = Sxx((i-1)*2+1)*sinb^2-2*Sxy((i-1)*2+1)*sinb*cosb + Syy((i-1)*2+1)*cosb^2;
    Tt((i-1)*2+1) = -(Sxx((i-1)*2+1) - Syy((i-1)*2+1))*sinb*cosb+ Sxy((i-1)*2+1)*(cosb^2 - sinb^2);
    Sn(i*2) = Sxx(i*2)*sinb^2-2*Sxy(i*2)*sinb*cosb + Syy(i*2)*cosb^2;
    Tt(i*2) = -(Sxx(i*2) - Syy(i*2))*sinb*cosb+ Sxy(i*2)*(cosb^2 - sinb^2);
end
for i = 1 : nEle
    nR = ConnList(i,4);
    nL = ConnList(i,3);
    if nR < 0.1
        beta1 = asind(AllEle(nL,5));
        beta2 = asind(AllEle(i,5));
        a2 = AllEle(i,7)/2;
        a1 = AllEle(nL,7)/2;
        ds = a2*cosd(beta2-beta1)+a1;
        
        us2 = Us((i-1)*2+1)*cosd(beta2 - beta1) - Un((i-1)*2+1)*sind(beta2-beta1);
        us1 = Us((nL-1)*2+1);
        duds = (us2 - us1)/ds;
        St((i-1)*2+1) = 2*Mat.G*1e3/(1-Mat.miu)*duds+Mat.miu/(1-Mat.miu)*Sn((i-1)*2+1);
        
        us2 = Us(i*2)*cosd(beta2 - beta1) - Un(i*2)*sind(beta2-beta1);
        us1 = Us(nL*2);
        duds = (us2 - us1)/ds;
        St(i*2) = 2*Mat.G*1e3/(1-Mat.miu)*duds+Mat.miu/(1-Mat.miu)*Sn(i*2);
    end
    if nL < 0.1
        beta1 = asind(AllEle(i,5));
        beta2 = asind(AllEle(nR,5));
        a1 = AllEle(i,7)/2;
        a2 = AllEle(nR,7)/2;
        ds = a2*cosd(beta2-beta1)+a1;
        
        us2 = Us((nR-1)*2+1)*cosd(beta2 - beta1) - Un((nR-1)*2+1)*sind(beta2-beta1);
        us1 = Us((i-1)*2+1);
        duds = (us2 - us1)/ds;
        St((i-1)*2+1) = 2*Mat.G*1e3/(1-Mat.miu)*duds+Mat.miu/(1-Mat.miu)*Sn((i-1)*2+1);
        
        us2 = Us(nR*2)*cosd(beta2 - beta1) - Un(nR*2)*sind(beta2-beta1);
        us1 = Us(i*2);
        duds = (us2 - us1)/ds;
        St(i*2) = 2*Mat.G*1e3/(1-Mat.miu)*duds+Mat.miu/(1-Mat.miu)*Sn(i*2);
    end
    if nL > 0.1 && nR > 0.1
        beta2 = asind(AllEle(i,5));
        beta3 = asind(AllEle(nR,5));
        beta1 = asind(AllEle(nL,5));
        a3 = AllEle(nR,7)/2;
        a2 = AllEle(i,7)/2;
        a1 = AllEle(nL,7)/2;
        
        ds = a1*cosd(beta2-beta1)+2*a2+a3*cosd(beta3-beta2);
        us2 = Us((nR-1)*2+1)*cosd(beta3 - beta2) - Us((nL-1)*2+1)*cosd(beta2-beta1);
        us1 = Un((nR-1)*2+1)*sind(beta3 - beta2) + Un((nL-1)*2+1)*sind(beta2-beta1);
        duds = (us2 - us1)/ds;
        St((i-1)*2+1) = 2*Mat.G*1e3/(1-Mat.miu)*duds+Mat.miu/(1-Mat.miu)*Sn((i-1)*2+1);
        
        us2 = Us(nR*2)*cosd(beta3 - beta2) - Us(nL*2)*cosd(beta2-beta1);
        us1 = Un(nR*2)*sind(beta3 - beta2) + Un(nL*2)*sind(beta2-beta1);
        duds = (us2 - us1)/ds;
        St(i*2) = 2*Mat.G*1e3/(1-Mat.miu)*duds+Mat.miu/(1-Mat.miu)*Sn(i*2);
    end
end


axis(PicScale);

amp = 1e-4/max(abs(Ux)+abs(Uy));
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
figure(5);
%amp = 0.05/max(abs(Sxx)+abs(Syy)+abs(Sxy));
hold off;
STR = St/5;
for i = 1 : nEle
    plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b','Linewidth',1.5);
    hold on;
    if STR((i-1)*2+1)> 0
        plot([xx((i-1)*2+1),xx((i-1)*2+1)-abs(STR((i-1)*2+1))*amp/2*AllEle(i,5)],[yy((i-1)*2+1),yy((i-1)*2+1)+abs(STR((i-1)*2+1))*amp/2*AllEle(i,6)],'k','LineWidth',1.5);
    else
        plot([xx((i-1)*2+1),xx((i-1)*2+1)-abs(STR((i-1)*2+1))*amp/2*AllEle(i,5)],[yy((i-1)*2+1),yy((i-1)*2+1)+abs(STR((i-1)*2+1))*amp/2*AllEle(i,6)],'r','LineWidth',1.5);
    end
    if STR(i*2)> 0
        plot([xx(i*2),xx(i*2)+abs(STR(i*2))*amp/2*AllEle(i,5)],[yy(i*2),yy(i*2)-abs(STR(i*2))*amp/2*AllEle(i,6)],'k--','LineWidth',1.5);
    else
        plot([xx(i*2),xx(i*2)+abs(STR(i*2))*amp/2*AllEle(i,5)],[yy(i*2),yy(i*2)-abs(STR(i*2))*amp/2*AllEle(i,6)],'r--','LineWidth',1.5);
    end
end
axis(PicScale);
axis equal;
title('Near Fracture Tangential Stress','Fontsize',13);
figure(16);
%amp = 0.05/max(abs(Sxx)+abs(Syy)+abs(Sxy));
hold off;
STR = Sn/5;
for i = 1 : nEle
    plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b','Linewidth',1.5);
    hold on;
    if STR((i-1)*2+1)> 0
        plot([xx((i-1)*2+1),xx((i-1)*2+1)-abs(STR((i-1)*2+1))*amp/2*AllEle(i,5)],[yy((i-1)*2+1),yy((i-1)*2+1)+abs(STR((i-1)*2+1))*amp/2*AllEle(i,6)],'k','LineWidth',1.5);
    else
        plot([xx((i-1)*2+1),xx((i-1)*2+1)-abs(STR((i-1)*2+1))*amp/2*AllEle(i,5)],[yy((i-1)*2+1),yy((i-1)*2+1)+abs(STR((i-1)*2+1))*amp/2*AllEle(i,6)],'r','LineWidth',1.5);
    end
    if STR(i*2)> 0
        plot([xx(i*2),xx(i*2)+abs(STR(i*2))*amp/2*AllEle(i,5)],[yy(i*2),yy(i*2)-abs(STR(i*2))*amp/2*AllEle(i,6)],'k--','LineWidth',1.5);
    else
        plot([xx(i*2),xx(i*2)+abs(STR(i*2))*amp/2*AllEle(i,5)],[yy(i*2),yy(i*2)-abs(STR(i*2))*amp/2*AllEle(i,6)],'r--','LineWidth',1.5);
    end
end
axis(PicScale);
axis equal;
title('Near Fracture Normal Stress','Fontsize',13);
figure(7);
%amp = 0.05/max(abs(Sxx)+abs(Syy)+abs(Sxy));
hold off;
STR = Tt/5;
for i = 1 : nEle
    plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b','Linewidth',1.5);
    hold on;
    if STR((i-1)*2+1)> 0
        plot([xx((i-1)*2+1),xx((i-1)*2+1)-abs(STR((i-1)*2+1))*amp/2*AllEle(i,5)],[yy((i-1)*2+1),yy((i-1)*2+1)+abs(STR((i-1)*2+1))*amp/2*AllEle(i,6)],'k','LineWidth',1.5);
    else
        plot([xx((i-1)*2+1),xx((i-1)*2+1)-abs(STR((i-1)*2+1))*amp/2*AllEle(i,5)],[yy((i-1)*2+1),yy((i-1)*2+1)+abs(STR((i-1)*2+1))*amp/2*AllEle(i,6)],'r','LineWidth',1.5);
    end
    if STR(i*2)> 0
        plot([xx(i*2),xx(i*2)+abs(STR(i*2))*amp/2*AllEle(i,5)],[yy(i*2),yy(i*2)-abs(STR(i*2))*amp/2*AllEle(i,6)],'k--','LineWidth',1.5);
    else
        plot([xx(i*2),xx(i*2)+abs(STR(i*2))*amp/2*AllEle(i,5)],[yy(i*2),yy(i*2)-abs(STR(i*2))*amp/2*AllEle(i,6)],'r--','LineWidth',1.5);
    end
end
axis(PicScale);
axis equal;
title('Near Fracture Shear Stress','Fontsize',13);
clear WatchLine;
nele = length(xx);
WatchLine(:,1) = xx(2:nele-1);
WatchLine(:,2) = yy(2:nele-1);
figure(8);
[UxL,UyL,SxxL,SyyL,SxyL,~,~] = CalcPointStress_H(Mat,DD,WatchLine(:,1),WatchLine(:,2),nFracture,Fractures,nEle,AllEle,ConnList);
SxxL = SxxL + SigmaX;
SyyL = SyyL + SigmaY;
SxyL = SxyL + SigmaXY;

[ne,~ ] = size(WatchLine);
seg = [WatchLine(1,1),WatchLine(1,2),WatchLine(2,1),WatchLine(2,2)];
dis = sqrt((seg(3)-seg(1))^2+(seg(4)-seg(2))^2);
cosb = (seg(3)-seg(1))/dis;
sinb = (seg(4)-seg(2))/dis;
UsL = zeros(ne,1);
UnL = zeros(ne,1);
SnL = zeros(ne,1);
StL = zeros(ne,1);
TtL = zeros(ne,1);

for i = 1 : ne
    UsL(i) = UxL(i) *cosb + UyL(i)*sinb;
    UnL(i) = -UxL(i) *sinb + UyL(i)*cosb;
    SnL(i) = SxxL(i)*sinb^2-2*SxyL(i)*sinb*cosb + SyyL(i)*cosb^2;
    TtL(i) = -(SxxL(i) - SyyL(i))*sinb*cosb+ SxyL(i)*(cosb^2 - sinb^2);
    StL(i) = SxxL(i)*cosb^2+2*SxyL(i)*sinb*cosb + SyyL(i)*sinb^2;
end
a = 0.5*(WatchLine(2,1)-WatchLine(1,1));
for i = 1 : ne
    ds = a*2;
    if i==1
        us2 = UsL(i+1) ;
        us1 = UsL(i);
        duds = (us2 - us1)/ds;
        StL(i) = 2*Mat.G*1e3/(1-Mat.miu)*duds+Mat.miu/(1-Mat.miu)*SnL(i);
    end
    if i == ne
        us2 = UsL(i) ;
        us1 = UsL(i-1);
        duds = (us2 - us1)/ds;
        StL(i) = 2*Mat.G*1e3/(1-Mat.miu)*duds+Mat.miu/(1-Mat.miu)*SnL(i);
    end
    if i < ne && i > 1
        ds = a*4;
        us2 = UsL(i+1)- UsL(i-1);
        us1 = 0;
        duds = (us2 - us1)/ds;
        StL(i) = 2*Mat.G*1e3/(1-Mat.miu)*duds+Mat.miu/(1-Mat.miu)*SnL(i);
    end
end
hold off;
xlim([-0.2 0.2])
%axis(PicScale);
plot(WatchLine(:,2),StL','.');
hold on;
title('Near NF Stress','Fontsize',13);
hold on;
plot(WatchLine(:,2),TtL,'k.');
plot(WatchLine(:,2),SnL,'g.');

legend('Tangential Stress','Shear Stress','Normal Stress');
for i = 1 : nwell
    plot([well{i}.heel(1),well{i}.toe(1)],[well{i}.heel(2),well{i}.toe(2)],'k','Linewidth',2);
    for j = 1 : well{i}.nPerf
        num = well{i}.Perfindex(j);
        plot([AllEle(num,1) AllEle(num,3)],[AllEle(num,2) AllEle(num,4)],'r','Linewidth',3);
        plot(well{i}.Perf(j,1),well{i}.Perf(j,2),'wo');
        plot(well{i}.Perf(j,1),well{i}.Perf(j,2),'g.','Markersize',20);
        hold on;
    end
end
title(['Stress Along Fractures @ Time ',num2str(curT)/60,'min'],'Fontsize',14);
saveas(16,[FILEPATH,num2str(nt),'CP.fig']);

end
