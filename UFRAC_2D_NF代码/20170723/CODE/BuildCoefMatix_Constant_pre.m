%***********************************************************
%********* Compute the Coefficient Matrix for  DD  *********
%*******            Huiying Tang @ 2015     check kink      *********
%**********************************************************
function [extraBC,CM] = BuildCoefMatix_Constant_pre(Mat,nAct,FracHeight,AllEle,TipStates)
% Higher order elements (Langrange 2nd order elements),u
% to note that nEle is a vector not a scalar
% Add Kink Element Consideration
global alpha;
global beta;
global UseHeight;
global UseTipE;

% Find Active Grids
extraBC = zeros(nAct*2,1);
G = Mat.G*1e9;
miu = Mat.miu;
GM = zeros(nAct,nAct);

if UseTipE > 0.1
    disp('Crack Tip Element is used');
else
    disp('Crack Tip Element is NOT used');%
end
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
        xjm = AllEle(j,8);
        yjm = AllEle(j,9);
        dij = sqrt((x-xjm)^2 + ( y - yjm)^2);
        if UseHeight > 0.1
            GM(i,j) = 1 - dij^beta/(dij^2 + (FracHeight/alpha)^2)^(beta/2);
        end
        isTip = 0;
        % Decide the tip element is along or opposite the direction of
        % fractures
        %ÖØÐÂÅÐ¶ÏisTip
        if abs(TipStates(j)) > 0.1
            if j > 1.1
                type = 1;
            else
                type = 2;
            end
            if type > 1.9
                isTip = 1;
            else
                isTip = 2;
            end
        end
        % ------------------------------ Not at Tip
        % ---------------------------------------------------
        % alpha is the angle of normal directions at ith element
        % For non-tip elements, Axn .. has three dimensions of Axn[3]
        % from three neightboured elements
        FracLoc = [AllEle(j,1),AllEle(j,2),AllEle(j,3),AllEle(j,4)];
        length = AllEle(j,7);
        if isTip < 0.9
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
        else
            % -------------Tip element
            % ---------------------------------
            FracLoc = [AllEle(j,1),AllEle(j,2),AllEle(j,3),AllEle(j,4)];
            lengthe = AllEle(j,7);
            if UseTipE  > 0.1
                [~,~,~,~,Axn,Axs,Ayn,Ays,Axyn,Axys]= CalcCoefficientMatrix_C(isTip,G,miu,FracLoc,x,y,AllEle(j,5),AllEle(j,6),lengthe);
                
            else
                [~,~,~,~,Axn,Axs,Ayn,Ays,Axyn,Axys]= CalcCoefficientMatrix_C(0,G,miu,FracLoc,x,y,AllEle(j,5),AllEle(j,6),lengthe);
            end
            temp1 = sinalp*cosalp*(Axs - Ays) - (cosalp^2 - sinalp^2)*Axys;
            temp2 = sinalp*cosalp*(Axn - Ayn) - (cosalp^2 - sinalp^2)*Axyn;
            CM(i,j) = CM(i,j)+temp1;
            CM(i,j+nAct) =CM(i,j+nAct) +temp2;
            % normal stress on element i
            temp1 = cosalp^2*Axs + 2*sinalp*cosalp*Axys + sinalp^2 * Ays;
            temp2 = cosalp^2*Axn + 2*sinalp*cosalp*Axyn + sinalp^2 * Ayn;
            CM(i+nAct,j) =CM(i+nAct,j)+ temp1;
            CM(i+nAct,j+nAct) = CM(i+nAct,j+nAct)+temp2;
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
