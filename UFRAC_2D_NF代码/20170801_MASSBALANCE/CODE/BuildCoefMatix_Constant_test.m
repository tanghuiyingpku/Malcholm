%***********************************************************
%********* Compute the Coefficient Matrix for  DD  *********
%*******     test for virtual element    *********
%**********************************************************
function [extraBC,CM] = BuildCoefMatix_Constant_test(nAct,Mat,AllEle)
% Higher order elements (Langrange 2nd order elements),u
% to note that nEle is a vector not a scalar
% Add Kink Element Consideration
global GM_global
global alpha;
global beta;
global UseHeight;
global CM_global
global MaxEle IndexInv;
FracHeight =Mat.h;
% Find Active Grids
extraBC = zeros(nAct*2,1);
G = Mat.G*1e9;
miu = Mat.miu;
GM = zeros(nAct,nAct);

if UseHeight > 0.1
    disp('Use Height Correction');
else
    disp('Not Use Height Correction');
end
CM = zeros(nAct*2, nAct*2);

% BC --[ shear stress...;normal stres....] DD-[shear opening,...; normal
% opening....]

% if the normal stress and shear stress are known
% Stress boundary condition

for i = 1 : nAct
    sinalp = AllEle(i,6);  %cosbet(i);
    cosalp = -AllEle(i,5); %-sinbet(i);
    x = (AllEle(i,1) + AllEle(i,3))/2;
    y = (AllEle(i,2) + AllEle(i,4))/2;
    %&& CM_marker(i)<2.9 menas its not a new element
    for j = 1 : nAct
        if i < nAct -0.1 && j < nAct-0.1
            CM(i,j) = CM_global(IndexInv(i),IndexInv(j));
            CM(i,j+nAct) = CM_global(IndexInv(i),IndexInv(j)+MaxEle);
            CM(i+nAct,j) = CM_global(IndexInv(i)+MaxEle,IndexInv(j));
            CM(i+nAct,j+nAct) = CM_global(IndexInv(i)+MaxEle,IndexInv(j)+MaxEle);
            if UseHeight > 0.1
                GM(i,j) = GM_global(IndexInv(i),IndexInv(j));
            end
            continue;
        else
            xjm = AllEle(j,8);
            yjm = AllEle(j,9);
            dij = sqrt((x-xjm)^2 + ( y - yjm)^2);
            if UseHeight > 0.1
                GM(i,j) = 1 - dij^beta/(dij^2 + (FracHeight/alpha)^2)^(beta/2);
            end
            % ------------------------------ Not at Tip
            % ---------------------------------------------------
            % alpha is the angle of normal directions at ith element
            % For non-tip elements, Axn .. has three dimensions of Axn[3]
            % from three neightboured elements
            FracLoc = [AllEle(j,1),AllEle(j,2),AllEle(j,3),AllEle(j,4)];
            length = AllEle(j,7);
            [~,~,~,~,Axn,Axs,Ayn,Ays,Axyn,Axys]= CalcCoefficientMatrix_C(0,G,miu,FracLoc,x,y,AllEle(j,5),AllEle(j,6),length);
            % transfer to the element i surface
            % shear stress on element i
            temp1 = sinalp*cosalp*(Axs - Ays) - (cosalp^2 - sinalp^2)*Axys;
            temp2 = sinalp*cosalp*(Axn - Ayn) - (cosalp^2 - sinalp^2)*Axyn;
            CM(i,j) = temp1;
            CM(i,j+nAct) = temp2;
            % normal stress on element i
            temp1 = cosalp^2*Axs + 2*sinalp*cosalp*Axys + sinalp^2 * Ays;
            temp2 = cosalp^2*Axn + 2*sinalp*cosalp*Axyn + sinalp^2 * Ayn;
            CM(i+nAct,j) = temp1;
            CM(i+nAct,j+nAct) = temp2;
        end
    end
end
if UseHeight > 0.1
    for i = 1 : nAct
        for j = 1 : nAct;
            CM(i,j)  = CM(i,j)*GM(i,j);
            CM(i,j+nAct)  = CM(i,j+nAct)*GM(i,j);
            CM(i+nAct,j)  = CM(i+nAct,j)*GM(i,j);
            CM(i+nAct,j+nAct)  = CM(i+nAct,j+nAct)*GM(i,j);
        end
    end
end
end
