function [Sxx,Syy,Sxy,X,Y] = CalcFieldStress_C(Mat,DD,xscale,yscale,AllEle,ConnList,SigmaX,SigmaY,SigmaXY)
% Calc Field Stress Surrounding Fractures
% xscale [x0,x1] the box size yscale

global UseTipE;
global UseHeight;
global alpha;
global beta;
global FracHeight;% unit:m
global IndexInv nAct;
global TipStates;
% Distance to Fracture surfaces
dis = 1e-3;
% Calculate Isolated range
if UseHeight > 0.1
    GM = zeros(nAct,1);
end
nDD = length(DD);
nEle = nDD/2;
% Write Fracture File
% Frac2 = AllEle(:,1:4);
% Frac2(:,1) = (Frac2(:,1) - xscale(1))/(xscale(2)-xscale(1));
% Frac2(:,3) = (Frac2(:,3) - xscale(1))/(xscale(2)-xscale(1));
% Frac2(:,2) = (Frac2(:,2) - yscale(1))/(yscale(2)-yscale(1));
% Frac2(:,4) = (Frac2(:,4) - yscale(1))/(yscale(2)-yscale(1));
% 
% refineN = 5;
% Frac = RefineFrac(nEle,Frac2,refineN);
% nEle2 = nEle*refineN;
% Frac = [Frac,ones(nEle2,1)];
% save Frac.dat -ascii Frac;
% foi = pwd;
% fracfile = [foi,'\Frac.dat'];
% % Trangle Center Points
% cd .././TriangleGridding/
% Point = GenerateTriangle(xscale(2)-xscale(1),yscale(2) - yscale(1),fracfile);
% cd ./CODE
% Point(:,1) = Point(:,1) + xscale(1);
% Point(:,2) = Point(:,2) + yscale(1);
N = 30;
xx = linspace(xscale(1),xscale(2),N);
yy = linspace(yscale(1),yscale(2),N);
for i = 1 : N
    Point((i-1)*N+1:i*N,:) = [xx',yy(i)*ones(N,1)];
end
[nGrid,~] = size(Point);
Sxx = zeros(nGrid, 1);
Syy = Sxx;
Sxy = Sxx;
X = Sxx;
Y = Sxx;
% to note that nEle is a vector not a scalar
G = Mat.G*1e3;
miu = Mat.miu;
% BC --[ shear stress...;normal stres....] DD-[shear opening,...; normal
% opening....]
tolDis = AllEle(1,7)/6;
tolDisTIP= AllEle(1,7)/3;
for i = 1 : nGrid
    x = Point(i,1);
    y = Point(i,2);
    Disc = 1e6;
    index = 0;
    for j = 1 : nEle
        dist = CalculateDis([x,y],AllEle(j,1:2));
        if dist < Disc
            Disc = dist;
            index = j;
        end
        dist = CalculateDis([x,y],AllEle(j,3:4));
        if dist < Disc
            Disc = dist;
            index = j;
        end
        if dist < Disc
            Disc = dist;
            index = j;
        end
        dist = CalculateDis([x,y],AllEle(j,8:9));
        if dist < Disc
            Disc = dist;
            index = j;
        end
    end
    if Disc < tolDis || (TipStates(index) > 0.1 && Disc < tolDisTIP)
        X(i) = x;
        Y(i) = y;
        Syy(i) = NaN;
        Sxy(i) = NaN;
        Sxx(i) = NaN;
        continue;
    end
    
    
    Rowy = zeros(nDD,1)';
    Rowxy = Rowy;
    Rowx = Rowy;
    for ne = 1 : nAct
        xjm = AllEle(IndexInv(ne),7);
        yjm = AllEle(IndexInv(ne),8);
        dij = sqrt((x-xjm)^2 + ( y - yjm)^2);
        if UseHeight > 0.1
            GM(ne) = 1 - dij^beta/(dij^2 + (FracHeight/alpha)^2)^(beta/2);
        end
        isTip = 0;
        if abs(TipStates(IndexInv(ne))) > 0.1
            [type,~] = TipType(IndexInv(ne),AllEle,ConnList);
            if type > 1.9
                isTip = 1;
            else
                isTip = 2;
            end
        end
        FracLoc = [AllEle(IndexInv(ne),1),AllEle(IndexInv(ne),2),AllEle(IndexInv(ne),3),AllEle(IndexInv(ne),4)];
        lengthe = AllEle(IndexInv(ne),7);
        sinbet = AllEle(IndexInv(ne),5);
        cosbet = AllEle(IndexInv(ne),6);
        if UseTipE < 0.1
            [~,~,~,~,Rowx(ne+nAct),Rowx(ne),Rowy(ne+nAct),Rowy(ne),Rowxy(ne+nAct),Rowxy(ne)]=CalcCoefficientMatrix_C(0,G,miu,FracLoc,x,y,sinbet,cosbet,lengthe);
        else
            [~,~,~,~,Rowx(ne+nAct),Rowx(ne),Rowy(ne+nAct),Rowy(ne),Rowxy(ne+nAct),Rowxy(ne)]=CalcCoefficientMatrix_C(isTip,G,miu,FracLoc,x,y,sinbet,cosbet,lengthe);
        end
    end
    if UseHeight > 0.1
        for nn = 1 : nAllEle
            Rowx(nn)  = Rowx(nn)*GM(nn);
            Rowx(nn+nAllEle)  = Rowx(nn+nAllEle)*GM(nn);
            Rowy(nn)  = Rowy(nn)*GM(nn);
            Rowy(nn+nAllEle)  = Rowy(nn+nAllEle)*GM(nn);
            Rowxy(nn)  = Rowxy(nn)*GM(nn);
            Rowxy(nn+nAllEle)  = Rowxy(nn+nAllEle)*GM(nn);
        end
    end
    X(i) = x;
    Y(i) = y;
    Syy(i) = Rowy * DD + SigmaY;
    Sxy(i) = Rowxy * DD + SigmaXY;
    Sxx(i) = Rowx *DD + SigmaX;
end
figure;
[XX,YY] = meshgrid(xx,yy);
Pm2 = griddata(X,Y,Sxx,XX,YY);
title('Sxx');
hold on;
surf(XX,YY,Pm2);
shading flat;
figure;
[XX,YY] = meshgrid(xx,yy);
Pm2 = griddata(X,Y,Sxy,XX,YY);
title('Sxy');
hold on;
surf(XX,YY,Pm2);
shading flat;
figure;
[XX,YY] = meshgrid(xx,yy);
Pm2 = griddata(X,Y,Syy,XX,YY);
% contourf(XX,YY,Pm2)
title('Syy');
hold on;
surf(XX,YY,Pm2);
shading flat;
[Sh,SH] = PrinStress(Sxx,Syy,Sxy);
figure
title('SH');
[XX,YY] = meshgrid(xx,yy);
Pm2 = griddata(X,Y,SH,XX,YY);
hold on;
surf(XX,YY,Pm2);
shading flat;
figure
title('Sh');
[XX,YY] = meshgrid(xx,yy);
Pm2 = griddata(X,Y,Sh,XX,YY);
hold on;
surf(XX,YY,Pm2);
shading flat;

end