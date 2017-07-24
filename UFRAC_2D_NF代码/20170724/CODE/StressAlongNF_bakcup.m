function  StressAlongNF_bakcup()
global  Index IndexInv nAct DD_global;
global  Mat  AllEle_global ConnList_global MaxEle nAllEle_global PresF_global;

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
eps = 1e-8;
% Draw NF
num = 0;

xx2 = zeros(500,1);
yy2 = zeros(500,1);
cosb = zeros(500,1);
sinb = zeros(500,1);
P = zeros(500,1);
for i = 1 : nAllEle_global
    if AllEle_global(i,10) > 1.1 
        num = num+ 1;
        xx2(num) = AllEle_global(i,8)-eps*AllEle_global(i,5);
        yy2(num) = AllEle_global(i,9)+eps*AllEle_global(i,6);
        cosb(num) =  AllEle_global(i,6);
        sinb(num) = AllEle_global(i,5);
        P(num) = PresF_global(i);
        num = num + 1;
        xx2(num) = AllEle_global(i,8)+eps*AllEle_global(i,5);
        yy2(num) = AllEle_global(i,9)-eps*AllEle_global(i,6);
        cosb(num) =  AllEle_global(i,6);
        sinb(num) = AllEle_global(i,5);
        P(num) = PresF_global(i);
    end
end
xx2 = xx2(1:num);
yy2 = yy2(1:num);
sigmaN = zeros(num,1);
tau = zeros(num,1);
St = zeros(num,1);
Us = zeros(num,1);
Un = zeros(num,1);
InsituS = zeros(num,1);

[Ux,Uy,Sx,Sy,Sxy,~,~] =  CalcPointStress_C_local(Mat,DD,xx2,yy2,AllEle);
for i = 1 : num
    Sxxi = Sx(i) + Mat.Sxx;
    Syyi = Sy(i) + Mat.Syy;
    Sxyi = Sxy(i) + Mat.Sxy;
    cosi = -sinb(i);
    sini = cosb(i);
    Us(i) = Ux(i) *cosi + Uy(i)*sini;
    Un(i) = -Ux(i) *sini + Uy(i)*cosi;
    InsituS(i) = (Mat.Sxx - Mat.Syy) * cosi * sini - Mat.Sxy*(cosi^2 - sini^2);
    sigmaN(i) = (Sxxi * cosi^2 + 2 * Sxyi * sini * cosi + Syyi * sini^2);
    tau(i) = (Sxxi- Syyi)*sini*cosi-Sxyi*(cosi^2 - sini^2);
end
% Tangential

figure(55);
hold off;

amp= 1E-2;
for i = 1 : num/2
    cosi = cosb(i*2);
    sini = sinb(i*2);
    plot(xx2(i*2),yy2(i*2),'ko');
    hold on;
    if sigmaN((i-1)*2+1)> 0
        plot([xx2((i-1)*2+1),xx2((i-1)*2+1)-abs(sigmaN((i-1)*2+1))*amp/2*sini],[yy2((i-1)*2+1),yy2((i-1)*2+1)+abs(sigmaN((i-1)*2+1))*amp/2*cosi],'k','LineWidth',1.5);
    else
        plot([xx2((i-1)*2+1),xx2((i-1)*2+1)-abs(sigmaN((i-1)*2+1))*amp/2*sini],[yy2((i-1)*2+1),yy2((i-1)*2+1)+abs(sigmaN((i-1)*2+1))*amp/2*cosi],'r','LineWidth',1.5);
    end
     plot([xx2((i-1)*2+1),xx2((i-1)*2+1)-abs(P((i-1)*2+1))*amp/2*sini],[yy2((i-1)*2+1),yy2((i-1)*2+1)+abs(P((i-1)*2+1))*amp/2*cosi],'b--','LineWidth',1.5);

    if sigmaN(i*2)> 0
        plot([xx2(i*2),xx2(i*2)+abs(sigmaN(i*2))*amp/2*sini],[yy2(i*2),yy2(i*2)-abs(sigmaN(i*2))*amp/2*cosi],'k--','LineWidth',1.5);
    else
        plot([xx2(i*2),xx2(i*2)+abs(sigmaN(i*2))*amp/2*sini],[yy2(i*2),yy2(i*2)-abs(sigmaN(i*2))*amp/2*cosi],'r--','LineWidth',1.5);
    end
end
title('Normal Stress');
axis equal;

figure(56);
hold off;
for i = 1 : num/2
    cosi = cosb(i*2);
    sini = sinb(i*2);
     plot(xx2(i*2),yy2(i*2),'ko');
     hold on;
    if tau((i-1)*2+1)> 0
        plot([xx2((i-1)*2+1),xx2((i-1)*2+1)-abs(tau((i-1)*2+1))*amp/2*sini],[yy2((i-1)*2+1),yy2((i-1)*2+1)+abs(tau((i-1)*2+1))*amp/2*cosi],'k','LineWidth',1.5);
    else
        plot([xx2((i-1)*2+1),xx2((i-1)*2+1)-abs(tau((i-1)*2+1))*amp/2*sini],[yy2((i-1)*2+1),yy2((i-1)*2+1)+abs(tau((i-1)*2+1))*amp/2*cosi],'r','LineWidth',1.5);
    end
    if tau(i*2)> 0
        plot([xx2(i*2),xx2(i*2)+abs(tau(i*2))*amp/2*sini],[yy2(i*2),yy2(i*2)-abs(tau(i*2))*amp/2*cosi],'k--','LineWidth',1.5);
    else
        plot([xx2(i*2),xx2(i*2)+abs(tau(i*2))*amp/2*sini],[yy2(i*2),yy2(i*2)-abs(tau(i*2))*amp/2*cosi],'r--','LineWidth',1.5);
    end
end
title('Shear Stress');
axis equal;
% 
% figure(57);
% for i = 1 : num/2
%     hold on;
%     cosi = cosb(i*2);
%     sini = sinb(i*2);
%     if Tn((i-1)*2+1)> 0
%         plot([xx2((i-1)*2+1),xx2((i-1)*2+1)-abs(Tn((i-1)*2+1))*amp/2*sini],[yy2((i-1)*2+1),yy2((i-1)*2+1)+abs(Tn((i-1)*2+1))*amp/2*cosi],'k','LineWidth',1.5);
%     else
%         plot([xx2((i-1)*2+1),xx2((i-1)*2+1)-abs(Tn((i-1)*2+1))*amp/2*sini],[yy2((i-1)*2+1),yy2((i-1)*2+1)+abs(Tn((i-1)*2+1))*amp/2*cosi],'r','LineWidth',1.5);
%     end
%     if Tn(i*2)> 0
%         plot([xx2(i*2),xx2(i*2)+abs(Tn(i*2))*amp/2*sini],[yy2(i*2),yy2(i*2)-abs(Tn(i*2))*amp/2*cosi],'k--','LineWidth',1.5);
%     else
%         plot([xx2(i*2),xx2(i*2)+abs(Tn(i*2))*amp/2*sini],[yy2(i*2),yy2(i*2)-abs(Tn(i*2))*amp/2*cosi],'r--','LineWidth',1.5);
%     end
% end
% title('Tangential Stress');

clear sigmaN
clear tau;
clear xx2;
clear yy2;

end

