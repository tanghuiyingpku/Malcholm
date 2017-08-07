function [Sxx,Syy,Sxy,X,Y] = CalcFieldStress_C_GLOBAL(nt)
% Calc Field Stress Surrounding Fractures
% xscale [x0,x1] the box size yscale
global ConnList_global AllEle_global Mat;
global DD_global MaxEle;
global UseTipE;
global UseHeight;
global alpha;
global beta;
global FILEPATH;% unit:m
global IndexInv nAct;
global TipStates PicScale;
Kundrain = Mat.E;
Kdrain = Mat.E*0.9;
biot = 1;
UseTipE = 0;
DD = zeros(nAct*2,1);
ConnList = zeros(nAct,8);
AllEle = AllEle_global;
for i = 1 : nAct
    DD(i) = DD_global(IndexInv(i));
    DD(i+nAct) = DD_global(IndexInv(i)+MaxEle);
  % AllEle(i,:) = AllEle_global(IndexInv(i),:);
    ConnList(i,:) = ConnList_global(IndexInv(i),:);
end
% Distance to Fracture surfaces
dis = 1e-3;
% Calculate Isolated range
if UseHeight > 0.1
    GM = zeros(nAct,1);
end
nDD = length(DD);
nEle = nDD/2;
SigmaX = Mat.Sxx ;
SigmaY = Mat.Syy ;
SigmaXY = Mat.Sxy ;
xscale(1) = PicScale(1);
xscale(2) = PicScale(2);
yscale(1) = PicScale(3);
yscale(2) = PicScale(4);

Nx = 80;
Ny = 100;
xx = linspace(xscale(1),xscale(2),Nx);
yy = linspace(yscale(1),yscale(2),Ny);
% for i = 1 : N
%     Point((i-1)*N+1:i*N,:) = [xx',yy(i)*ones(N,1)];
% end
% [nGrid,~] = size(Point);
nGrid = Nx*Ny;

Sxx = zeros(Nx, Ny);
Syy = Sxx;
Sxy = Sxx;
InduceP = Sxx;
X = Sxx;
Y = Sxx;
SH = Sxx;
Sh = Sxx;
% to note that nEle is a vector not a scalar
G = Mat.G*1e3;
miu = Mat.miu;
% BC --[ shear stress...;normal stres....] DD-[shear opening,...; normal
% opening....]
tolDis = AllEle(1,7)/6;
tolDisTIP= AllEle(1,7)/3;
for ii =1 : Nx
    for jj = 1 : Ny
        x = xx(ii);
        y = yy(jj);
        % for i = 1 : nGrid
        %     x = Point(i,1);
        %     y = Point(i,2);
       
        %     if Disc < tolDis || (TipStates(index) > 0.1 && Disc < tolDisTIP)
        %         X(i) = x;
        %         Y(i) = y;
        %         Syy(i) = NaN;
        %         Sxy(i) = NaN;
        %         Sxx(i) = NaN;
        %         continue;
        %     end
        Rowy = zeros(nDD,1)';
        Rowxy = Rowy;
        Rowx = Rowy;
        for ne = 1 : nAct
            xjm = AllEle(IndexInv(ne),7);
            yjm = AllEle(IndexInv(ne),8);
            dij = sqrt((x-xjm)^2 + ( y - yjm)^2);
            if UseHeight > 0.1
                GM(ne) = 1 - dij^beta/(dij^2 + (Mat.h/alpha)^2)^(beta/2);
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
            %         for nn = 1 : nAct
            %             Rowx(nn)  = Rowx(nn)*GM(nn);
            %             Rowx(nn+nAct)  = Rowx(nn+nAct)*GM(nn);
            %             Rowy(nn)  = Rowy(nn)*GM(nn);
            %             Rowy(nn+nAct)  = Rowy(nn+nAct)*GM(nn);
            %             Rowxy(nn)  = Rowxy(nn)*GM(nn);
            %             Rowxy(nn+nAct)  = Rowxy(nn+nAct)*GM(nn);
            %         end
        end
        X(ii,jj) = x;
        Y(ii,jj) = y;
        Syy(ii,jj) = Rowy * DD ;
        Sxy(ii,jj) = Rowxy * DD ;
        Sxx(ii,jj) = Rowx *DD ;
        %InduceP(ii,jj) = (Sxx(i) + Syy(i) + Mat.miu*(Sxx(i) + Syy(i)))*(1/3/Kundrain - 1/3/Kdrain)*Kdrain/biot;
        Syy(ii,jj) = Syy(ii,jj) + SigmaY;
        Sxy(ii,jj) =Sxy(ii,jj) + SigmaXY;
        Sxx(ii,jj) = Sxx(ii,jj)  + SigmaX;
        [aa,bb] = PrinStress(Sxx(ii,jj),Syy(ii,jj),Sxy(ii,jj));
        SH(ii,jj) = bb;
        Sh(ii,jj) = aa;
    end
end
figure(21);
title('Sxx');
hold on;
surf(X,Y,Sxx);
shading flat;
saveas(21,[FILEPATH,'Sxx.fig']);

figure(22);
title('Sxy');
hold on;
surf(X,Y,Sxy);
shading flat;
saveas(22,[FILEPATH,'SxY.fig']);

figure(23);
title('Syy');
hold on;
surf(X,Y,Syy);
shading flat;
saveas(23,[FILEPATH,'sYY.fig']);


figure(24);
title('SH');
hold on;
surf(X,Y,SH);
shading flat;
saveas(24,[FILEPATH,'ShMAX.fig']);

figure(25);
title('Sh');
hold on;
surf(X,Y,Sh);
shading flat;
saveas(25,[FILEPATH,'SHMIN.fig']);

end