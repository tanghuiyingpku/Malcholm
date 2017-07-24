function  StressAlongNF_tangential()
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
IndexNF = zeros(nAllEle_global,1);
for i = 1 : nAllEle_global
    if AllEle_global(i,10) > 1.1 
        num = num+ 1;
        xx2(num) = AllEle_global(i,8)-eps*AllEle_global(i,5);
        yy2(num) = AllEle_global(i,9)+eps*AllEle_global(i,6);
        cosb(num) =  AllEle_global(i,6);
        sinb(num) = AllEle_global(i,5);
        P(num) = PresF_global(i);
        num = num + 1;
        IndexNF(i) = num/2;
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
    Us(i) = Ux(i) *(sini) + Uy(i)*(-cosi);
    Un(i) = -Ux(i) *(-cosi) + Uy(i)*sini;
    InsituS(i) = (Mat.Sxx - Mat.Syy) * cosi * sini - Mat.Sxy*(cosi^2 - sini^2);
    sigmaN(i) = (Sxxi * cosi^2 + 2 * Sxyi * sini * cosi + Syyi * sini^2);
    tau(i) = (Sxxi- Syyi)*sini*cosi-Sxyi*(cosi^2 - sini^2);
end
% Tangential
num = 0;
for i = 1 : nAllEle_global
    if AllEle_global(i,10) > 1.1
        
        nR = ConnList_global(i,4);
        nL = ConnList_global(i,3);
        if nR < 0.1
            beta1 = asind(AllEle_global(nL,5));
            beta2 = asind(AllEle_global(i,5));
            a2 = AllEle_global(i,7)/2;
            a1 = AllEle_global(nL,7)/2;
            ds = a2*cosd(beta2-beta1)+a1;
            num = (IndexNF(i)-1)*2+1;
            numL = (IndexNF(nL)-1)*2+1;
            us2 = Us(num)*cosd(beta2 - beta1) - Un(num)*sind(beta2-beta1);
            us1 = Us(numL);
            duds = (us2 - us1)/ds;
            St(num) = 2*Mat.G*1e3/(1-Mat.miu)*duds+Mat.miu/(1-Mat.miu)*sigmaN(num);
            
            num = (IndexNF(i)-1)*2+2;
            numL = (IndexNF(nL)-1)*2+2;
            us2 = Us(num)*cosd(beta2 - beta1) - Un(num)*sind(beta2-beta1);
            us1 = Us(numL);
            duds = (us2 - us1)/ds;
            St(num) = 2*Mat.G*1e3/(1-Mat.miu)*duds+Mat.miu/(1-Mat.miu)*sigmaN(num);
        end
        if nL < 0.1
            beta1 = asind(AllEle_global(i,5));
            beta2 = asind(AllEle_global(nR,5));
            a1 = AllEle_global(i,7)/2;
            a2 = AllEle_global(nR,7)/2;
            ds = a2*cosd(beta2-beta1)+a1;
            
            num = (IndexNF(i)-1)*2+1;
            numR = (IndexNF(nR)-1)*2+1;
            us2 = Us(numR)*cosd(beta2 - beta1) - Un(numR)*sind(beta2-beta1);
            us1 = Us(num);
            duds = (us2 - us1)/ds;
            St(num) = 2*Mat.G*1e3/(1-Mat.miu)*duds+Mat.miu/(1-Mat.miu)*sigmaN(num);
            
            num = (IndexNF(i)-1)*2+2;
            numR = (IndexNF(nR)-1)*2+2;
            us2 = Us(numR)*cosd(beta2 - beta1) - Un(numR)*sind(beta2-beta1);
            us1 = Us(num);
            duds = (us2 - us1)/ds;
            St(num) = 2*Mat.G*1e3/(1-Mat.miu)*duds+Mat.miu/(1-Mat.miu)*sigmaN(num);
        end
        if nL > 0.1 && nR > 0.1
            beta2 = asind(AllEle_global(i,5));
            beta3 = asind(AllEle_global(nR,5));
            beta1 = asind(AllEle_global(nL,5));
            a3 = AllEle_global(nR,7)/2;
            a2 = AllEle_global(i,7)/2;
            a1 = AllEle_global(nL,7)/2;
            
            numL = (IndexNF(nL)-1)*2+1;
            num = (IndexNF(i)-1)*2+1;
            numR = (IndexNF(nR)-1)*2+1;
            
            ds = a1*cosd(beta2-beta1)+2*a2+a3*cosd(beta3-beta2);
            us2 = Us(numR)*cosd(beta3 - beta2) - Us(numL)*cosd(beta2-beta1);
            us1 = Un(numR)*sind(beta3 - beta2) + Un(numL)*sind(beta2-beta1);
            duds = (us2 - us1)/ds;
            St(num) = 2*Mat.G*1e3/(1-Mat.miu)*duds+Mat.miu/(1-Mat.miu)*sigmaN(num);
            
            numL = (IndexNF(nL)-1)*2+2;
            num = (IndexNF(i)-1)*2+2;
            numR = (IndexNF(nR)-1)*2+2;
            
            us2 = Us(numR)*cosd(beta3 - beta2) - Us(numL)*cosd(beta2-beta1);
            us1 = Un(numR)*sind(beta3 - beta2) + Un(numL)*sind(beta2-beta1);
            duds = (us2 - us1)/ds;
            St(num) = 2*Mat.G*1e3/(1-Mat.miu)*duds+Mat.miu/(1-Mat.miu)*sigmaN(num);
        end
    end
end


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

St = Sx;
figure(57);
for i = 1 : num/2
    hold on;
    cosi = cosb(i*2);
    sini = sinb(i*2);
    if St((i-1)*2+1)> 0
        plot([xx2((i-1)*2+1),xx2((i-1)*2+1)-abs(St((i-1)*2+1))*amp/2*sini],[yy2((i-1)*2+1),yy2((i-1)*2+1)+abs(St((i-1)*2+1))*amp/2*cosi],'k','LineWidth',1.5);
    else
        plot([xx2((i-1)*2+1),xx2((i-1)*2+1)-abs(St((i-1)*2+1))*amp/2*sini],[yy2((i-1)*2+1),yy2((i-1)*2+1)+abs(St((i-1)*2+1))*amp/2*cosi],'r','LineWidth',1.5);
    end
    if St(i*2)> 0
        plot([xx2(i*2),xx2(i*2)+abs(St(i*2))*amp/2*sini],[yy2(i*2),yy2(i*2)-abs(St(i*2))*amp/2*cosi],'k--','LineWidth',1.5);
    else
        plot([xx2(i*2),xx2(i*2)+abs(St(i*2))*amp/2*sini],[yy2(i*2),yy2(i*2)-abs(St(i*2))*amp/2*cosi],'r--','LineWidth',1.5);
    end
end
title('Tangential Stress');
axis equal;

clear sigmaN
clear tau;
clear xx2;
clear yy2 cosb  sinb  P  IndexNF ;
clear  St  Us  InsituS Un

end

