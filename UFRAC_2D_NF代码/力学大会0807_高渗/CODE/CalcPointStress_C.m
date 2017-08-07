function [Ux,Uy,Sxx,Syy,Sxy,X,Y] = CalcPointStress_C(Mat,DD,xx,yy,AllEle,ConnList)
global UseTipE;
global UseHeight;
global alpha;
global beta; 
global IndexInv ;
global TipStates;
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
nGridy = length(yy);
mark = zeros(nGridx,nGridy);
% for i = 1 : nGridx
%     for j = 1 : nGridy
%         dot = [xx(i),yy(j)];
%         for k =1 : nDD/2
%             test = isDotIn(dot,AllEle((k),1:4));
%             if test > 0.1
%                 mark(i,j) =1;
%             end
%         end
%     end
% end
        
        
Sxx = zeros(nGridx, nGridy);
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
for i = 1 : nGridy
    y = yy(i);
    for j = 1 : nGridx
        x = xx(j);
       
        Rowy = zeros(nDD,1)';
        Rowxy = Rowy;
        Rowx = Rowy;
        Rowxu  = Rowy;
        Rowyu = Rowy;
        
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