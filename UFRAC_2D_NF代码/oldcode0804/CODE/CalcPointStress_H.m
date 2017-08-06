function [Ux,Uy,Sxx,Syy,Sxy,X,Y] = CalcPointStress_H(Mat,DD,xx,yy,nFracture,Fractures,nAllEle,AllEle,ConnList)
global UseTipE;
global TipStates IndexInv;
global UseHeight;
global alpha;
global beta;
global FracHeight;% unit:m
disp('      ');
disp('Calculate Points Stresses with High Order');
disp('Note!: The Stress is only the INDUCED STRESS!!');
disp('Real Stress shoule + In-situ Stress');
disp('      ');

if UseHeight > 0.1
    GM = zeros(nAllEle,1);
    disp('Use Height Correction');
else
    disp('Not Use Height Correction');
end
if UseTipE > 0.1
    disp('Crack Tip Element is used');
else
    disp('Crack Tip Element is NOT used');
end
nDD = length(DD);
nX = length(xx);
nY = length(yy);
Sxx = zeros(nX, 1);
Syy = Sxx;
Sxy = Sxx;
Ux = Sxx;
Uy = Sxx;
X = Sxx;
Y = Sxx;

% to note that nEle is a vector not a scalar
G = Mat.G*1e3;
miu = Mat.miu;
% BC --[ shear stress...;normal stres....] DD-[shear opening,...; normal
% opening....]
for i = 1 : nX
    x = xx(i);
    
    %for j = 1 : nY
    y = yy(i);
    Rowy = zeros(nDD,1)';
    Rowxy = Rowy;
    Rowx = Rowy;
    Rowxu  = Rowy;
    Rowyu = Rowy;
    
    for ne = 1 : nAllEle
        xjm = AllEle(ne,8);
        yjm = AllEle(ne,9);
        dij = sqrt((x-xjm)^2 + ( y - yjm)^2);
        if UseHeight > 0.1
            GM(ne) = 1 - dij^beta/(dij^2 + (FracHeight/alpha)^2)^(beta/2);
        end
        isTip = 0;
        if abs(TipStates(IndexInv(ne))) > 0.1
            [type,~] = TipType2(ne,AllEle,ConnList);
            if type > 1.9
                isTip = 1;
            else
                isTip = 2;
            end
        end
        if isTip < 0.9
            LeftEle = ConnList(ne,3);
            RightEle = ConnList(ne,4);
            if LeftEle < 0 || RightEle < 0
                keyboard;
            end
            Length = [AllEle(LeftEle,7),AllEle(ne,7),AllEle(RightEle,7)];
            CurEle =[AllEle(LeftEle,1:4);AllEle(ne,1:4);AllEle(RightEle,1:4)];
            [Auxn,Auxs,Auyn,Auys,Axn,Axs,Ayn,Ays,Axyn,Axys] = CalcCoefficientMatrix_H(isTip,CurEle,Length,AllEle(ne,5),AllEle(ne,6),G,miu,x,y);
            % transfer to the element i surface
            % xx  stress on element i
            Rowxu(ne) = Rowxu(ne) + Auxs(2);
            Rowxu(LeftEle) = Rowxu(LeftEle) + Auxs(1);
            Rowxu(RightEle) = Rowxu(RightEle) + Auxs(3);
            Rowxu(LeftEle+nAllEle) = Rowxu(LeftEle+nAllEle) + Auxn(1);
            Rowxu(ne+nAllEle) = Rowxu(ne+nAllEle) + Auxn(2);
            Rowxu(RightEle+nAllEle) = Rowxu(RightEle+nAllEle) + Auxn(3);
            
            Rowyu(LeftEle) = Rowyu(LeftEle) + Auys(1);
            Rowyu(ne) = Rowyu(ne) + Auys(2);
            Rowyu(RightEle) = Rowyu(RightEle) + Auys(3);
            Rowyu(LeftEle+nAllEle) = Rowyu(LeftEle+nAllEle) + Auyn(1);
            Rowyu(ne+nAllEle) = Rowyu(ne+nAllEle) + Auyn(2);
            Rowyu(RightEle+nAllEle) = Rowyu(RightEle+nAllEle) + Auyn(3);
            
            
            Rowx(LeftEle) = Rowx(LeftEle)+Axs(1);
            Rowx(ne) = Rowx(ne)+Axs(2);
            Rowx(RightEle) =Rowx(RightEle)+ Axs(3);
            Rowx(LeftEle+nAllEle) = Rowx(LeftEle+nAllEle)+Axn(1);
            Rowx(ne+nAllEle) =Rowx(ne+nAllEle) + Axn(2);
            Rowx(RightEle+nAllEle) =Rowx(RightEle+nAllEle)+ Axn(3);
            % yy stress on element i
            Rowy(LeftEle) = Rowy(LeftEle)+Ays(1);
            Rowy(ne) =Rowy(ne)+ Ays(2);
            Rowy(RightEle) =Rowy(RightEle)+ Ays(3);
            Rowy(LeftEle+nAllEle) = Rowy(LeftEle+nAllEle)+Ayn(1);
            Rowy(ne+nAllEle) = Rowy(ne+nAllEle)+Ayn(2);
            Rowy(RightEle+nAllEle) =Rowy(RightEle+nAllEle)+ Ayn(3);
            % Shear xy stresses on element i
            Rowxy(LeftEle) = Rowxy(LeftEle)+Axys(1);
            Rowxy(ne) =Rowxy(ne)+ Axys(2);
            Rowxy(RightEle) =Rowxy(RightEle)+ Axys(3);
            Rowxy(LeftEle+nAllEle) = Rowxy(LeftEle+nAllEle)+Axyn(1);
            Rowxy(ne+nAllEle) = Rowxy(ne+nAllEle)+Axyn(2);
            Rowxy(RightEle+nAllEle) =Rowxy(RightEle+nAllEle)+ Axyn(3);
        else
            FracLoc = [AllEle(ne,1),AllEle(ne,2),AllEle(ne,3),AllEle(ne,4)];
            lengthe = AllEle(ne,7);
            if UseTipE > 0.1
                [Auxn,Auxs,Auyn,Auys,Axn,Axs,Ayn,Ays,Axyn,Axys]=CalcCoefficientMatrix_C(isTip,G,miu,FracLoc,x,y,AllEle(ne,5),AllEle(ne,6),lengthe);
            else
                [Auxn,Auxs,Auyn,Auys,Axn,Axs,Ayn,Ays,Axyn,Axys]=CalcCoefficientMatrix_C(0,G,miu,FracLoc,x,y,AllEle(ne,5),AllEle(ne,6),lengthe);
            end
            
            Rowyu(ne+nAllEle) = Rowyu(ne+nAllEle) + Auyn;
            Rowyu(ne) = Rowyu(ne) + Auys;
            Rowxu(ne+nAllEle) = Rowxu(ne+nAllEle) + Auxn;
            Rowxu(ne) = Rowxu(ne) + Auxs;
            % xx  stress on element i
            Rowx(ne) = Rowx(ne)+Axs;
            Rowx(ne+nAllEle) =Rowx(ne+nAllEle) + Axn;
            % yy stress on element i
            Rowy(ne) =Rowy(ne)+ Ays;
            Rowy(ne+nAllEle) = Rowy(ne+nAllEle)+Ayn;
            % Shear xy stresses on element i
            Rowxy(ne) =Rowxy(ne)+ Axys;
            Rowxy(ne+nAllEle) = Rowxy(ne+nAllEle)+Axyn;
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
    Ux(i) = Rowxu*DD;
    Uy(i) = Rowyu*DD;
    Syy(i) = Rowy * DD;
    Sxy(i) = Rowxy * DD;
    Sxx(i) = Rowx *DD;
    %end
end
end





