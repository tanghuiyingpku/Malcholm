function [Ux,Uy,Sxx,Syy,Sxy,X,Y] = CalcPointStress_C_local_2(isActive,Mat,DD,xx,yy,AllEle)
global UseTipE;
global UseHeight;
global alpha;
global beta; 
% 
% disp('      ');
% disp('Calculate Points Stresses with Constant Element');
% disp('Note!: The Stress is only the INDUCED STRESS!!');
% disp('Real Stress shoule + In-situ Stress');
% disp('      ');
nDD = length(DD);
nAct = nDD/2;
if UseHeight > 0.1
    GM = zeros(nAct,1);
end

nGridx = length(xx);
mark = zeros(nGridx,1);
        
Sxx = zeros(nGridx, 1);
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
for i = 1 : nGridx
    for j = 1 : 1
        x = xx(i);
        y = yy(i);
        Rowy = zeros(nDD,1)';
        Rowxy = Rowy;
        Rowx = Rowy;
        Rowxu  = Rowy;
        Rowyu = Rowy;
        
        for ne = 1 : nAct
            xjm = AllEle(ne,7);
            yjm = AllEle(ne,8);
            dij = sqrt((x-xjm)^2 + ( y - yjm)^2);
            if UseHeight > 0.1
                GM(ne) = 1 - dij^beta/(dij^2 + (Mat.h/alpha)^2)^(beta/2);
            end
            isTip = 0;
            FracLoc = [AllEle(ne,1),AllEle(ne,2),AllEle(ne,3),AllEle(ne,4)];
            lengthe = AllEle(ne,7);
            sinbet = AllEle(ne,5);
            cosbet = AllEle(ne,6);
            if UseTipE < 0.1
                [Rowxu(ne+nAct),Rowxu(ne),Rowyu(ne+nAct),Rowyu(ne),Rowx(ne+nAct),Rowx(ne),Rowy(ne+nAct),Rowy(ne),Rowxy(ne+nAct),Rowxy(ne)]=CalcCoefficientMatrix_C(0,G,miu,FracLoc,x,y,sinbet,cosbet,lengthe);
            else
                [Rowxu(ne+nAct),Rowxu(ne),Rowyu(ne+nAct),Rowyu(ne),Rowx(ne+nAct),Rowx(ne),Rowy(ne+nAct),Rowy(ne),Rowxy(ne+nAct),Rowxy(ne)]=CalcCoefficientMatrix_C(isTip,G,miu,FracLoc,x,y,sinbet,cosbet,lengthe);
            end
        end
        if UseHeight > 0.1
            for nn = 1 : nAct
                Rowx(nn)  = Rowx(nn)*GM(nn);
                Rowx(nn+nAct)  = Rowx(nn+nAct)*GM(nn);
                Rowy(nn)  = Rowy(nn)*GM(nn);
                Rowy(nn+nAct)  = Rowy(nn+nAct)*GM(nn);
                Rowxy(nn)  = Rowxy(nn)*GM(nn);
                Rowxy(nn+nAct)  = Rowxy(nn+nAct)*GM(nn);
            end
        end
        X(i,j) = xx(j);
        Y(i,j) = yy(i);
        Ux(i,j) = Rowxu*DD;
        Uy(i,j) = Rowyu*DD;
        Syy(i,j) = Rowy * DD;
        Sxy(i,j) = Rowxy * DD;
        Sxx(i,j) = Rowx * DD;
        if mark(i,j) > 0.1
            Syy(i,j) = 0;
            Sxy(i,j) = 0;
            Sxx(i,j) = 0;
        end
    end
end
end