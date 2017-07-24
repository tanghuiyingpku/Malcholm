%***********************************************************
%********* Compute the Coefficient Matrix for  DD  *********
%*******            Huiying Tang @ 2015     check kink      *********
%**********************************************************
function [extraBC,CM] = BuildCoefMatix_Constant2(Mat,Fractures,AllEle,ConnList,isInit)
% Higher order elements (Langrange 2nd order elements),u
% to note that nEle is a vector not a scalar
% Add Kink Element Consideration
global GM_global
global alpha;
global beta;
global UseHeight;
global UseTipE;
global TipStates;
global CM_global
global MaxEle;
global IndexInv nAct nActOld
FracHeight = Mat.h;
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
    sinalp = AllEle(IndexInv(i),6);  %cosbet(i);
    cosalp = -AllEle(IndexInv(i),5); %-sinbet(i);
    x = (AllEle(IndexInv(i),1) + AllEle(IndexInv(i),3))/2;
    y = (AllEle(IndexInv(i),2) + AllEle(IndexInv(i),4))/2;
    %&& CM_marker(i)<2.9 menas its not a new element
    for j = 1 : nAct
        if (isInit < 0.1 )|| (i <= nActOld && j<= nActOld)
            CM(i,j) = CM_global(IndexInv(i),IndexInv(j));
            CM(i,j+nAct) = CM_global(IndexInv(i),IndexInv(j)+MaxEle);
            CM(i+nAct,j) = CM_global(IndexInv(i)+MaxEle,IndexInv(j));
            CM(i+nAct,j+nAct) = CM_global(IndexInv(i)+MaxEle,IndexInv(j)+MaxEle);
            if UseHeight > 0.1
                GM(i,j) = GM_global(IndexInv(i),IndexInv(j));
            end
            continue;
        end
        xjm = AllEle(IndexInv(j),8);
        yjm = AllEle(IndexInv(j),9);
        dij = sqrt((x-xjm)^2 + ( y - yjm)^2);
        if UseHeight > 0.1
            GM(i,j) = 1 - dij^beta/(dij^2 + (FracHeight/alpha)^2)^(beta/2);
        end
        isTip = 0;
        % Decide the tip element is along or opposite the direction of
        % fractures
        %重新判断isTip
        if abs(TipStates(IndexInv(j))) > 0.1
            [type,~] = TipType(IndexInv(j),AllEle,ConnList);
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
        FracLoc = [AllEle(IndexInv(j),1),AllEle(IndexInv(j),2),AllEle(IndexInv(j),3),AllEle(IndexInv(j),4)];
        length = AllEle(IndexInv(j),7);
        if isTip < 0.9
            [~,~,~,~,Axn,Axs,Ayn,Ays,Axyn,Axys]= CalcCoefficientMatrix_C(0,G,miu,FracLoc,x,y,AllEle(IndexInv(j),5),AllEle(IndexInv(j),6),length);
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
            FracLoc = [AllEle(IndexInv(j),1),AllEle(IndexInv(j),2),AllEle(IndexInv(j),3),AllEle(IndexInv(j),4)];
            lengthe = AllEle(IndexInv(j),7);
            if UseTipE  > 0.1
                [~,~,~,~,Axn,Axs,Ayn,Ays,Axyn,Axys]= CalcCoefficientMatrix_C(isTip,G,miu,FracLoc,x,y,AllEle(IndexInv(j),5),AllEle(IndexInv(j),6),lengthe);
                
            else
                [~,~,~,~,Axn,Axs,Ayn,Ays,Axyn,Axys]= CalcCoefficientMatrix_C(0,G,miu,FracLoc,x,y,AllEle(IndexInv(j),5),AllEle(IndexInv(j),6),lengthe);
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

for i = 1 : nAct
    iparent = ConnList(IndexInv(i),1);
    %%%%%%%%%%%%%%%%%Stick Model%%%%%%%%%%%%%%
    if abs(AllEle(IndexInv(i),10) - 2) < 1e-6
        Kn = Fractures{iparent}.Kn;
        Ks = Fractures{iparent}.Ks;
        CM(i,i) = CM(i,i) + Ks*1e6;
        CM(i+nAct,i+nAct) = CM(i+nAct,i+nAct) + Kn*1e6;
    end
    %%%%%%%%%%%%%%%%% Slip Model %%%%%%%%%%%%%%%%
    if abs(AllEle(IndexInv(i),10) - 4 ) < 1e-6
        Kn = Fractures{iparent}.Kn;
        coh = Fractures{iparent}.Coh;
        frc = Fractures{iparent}.Fang;
        CM(i+nAct,i+nAct) = CM(i+nAct,i+nAct) + Kn*1e6;
        %上一个时间步符号
        signS = sign(AllEle(IndexInv(i),11));
        extraBC(i) = extraBC(i) + signS*coh*1e6;
        CM(i,i+nAct) = CM(i,i+nAct) - signS*Kn*1e6*tand(frc);
    end
    for j = 1 : nAct
        if UseHeight > 0.1
            CM(i,j)  = CM(i,j)*GM(i,j);
            CM(i,j+nAct)  = CM(i,j+nAct)*GM(i,j);
            CM(i+nAct,j)  = CM(i+nAct,j)*GM(i,j);
            CM(i+nAct,j+nAct)  = CM(i+nAct,j+nAct)*GM(i,j);
        end
        if isInit > 0.1
            %Update Global Value
            if UseHeight > 0.1
                GM_global(IndexInv(i),IndexInv(j)) = GM(i,j);
            end
            CM_global(IndexInv(i),IndexInv(j)) = CM(i,j);
            CM_global(IndexInv(i),MaxEle+IndexInv(j)) = CM(i,j+nAct);
            CM_global(IndexInv(i)+MaxEle,IndexInv(j)) = CM(i+nAct,j);
            CM_global(IndexInv(i)+MaxEle,IndexInv(j)+MaxEle) = CM(i+nAct,j+nAct);
        end
    end    
end

end
