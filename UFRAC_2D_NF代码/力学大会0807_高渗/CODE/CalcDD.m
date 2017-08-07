function CalcDD()
global  IndexInv nAct nAllEle_global EleType;
global   Mat DD_global MaxEle;
global  PresF_global AllEle_global ConnList_global nFracture Fractures;
Pres = zeros(nAct,1);
for i = 1 : nAct
    Pres(i) = PresF_global(IndexInv(i))*1e6;
end
% Project the in-situ stresses to fractures;
Sxx = Mat.Sxx*1e6;
Syy = Mat.Syy*1e6;
Sxy = Mat.Sxy*1e6;

InsituS = zeros(nAct*2,1);
AllEle = zeros(nAct,11);
ConnList = zeros(nAct,8);
BC = zeros(nAct*2,1);
for i = 1 : nAct
    AllEle(i,:) = AllEle_global(IndexInv(i),:);
    ConnList(i,:) = ConnList_global(IndexInv(i),:);
    %alpha = theta + 90;
    sinalp = AllEle(i,6);  %cosbet(i);
    cosalp = -AllEle(i,5); %-sinbet(i);
    InsituS(i) = (Sxx - Syy) * cosalp * sinalp - Sxy*(cosalp^2 - sinalp^2);
    InsituS(i + nAct) = Sxx * cosalp^2 + Syy * sinalp^2+2*Sxy*sinalp*cosalp;
end
BC_insitu = InsituS;

for ii = 1 : nAct
    BC(ii) = BC(ii)-BC_insitu(ii);
    BC(nAct+ii) =BC(nAct+ii) -BC_insitu(nAct+ii) - Pres(ii);
end

if abs(EleType - 1) < 1E-6
    disp('Building CCD Influence Coefficient Matrix');
    disp('           ');
    [extraBC,CM]  = BuildCoefMatix_Constant(Mat,Fractures,AllEle_global(1:nAllEle_global,:),ConnList_global(1:nAllEle_global,:),1);
    %Unkown Displacements
    %Solve Linear System for Displacements
    disp('           ');
    disp('Matrix Solving...');
    disp('           ');
    BC = BC + extraBC;
    DD = CM\BC;
    %Analytical
    c = 0.02;
    E = Mat.E;
    dp = Pres(1) + Mat.Sxx;
    dp = dp/1e6;
    w  = 4*dp/E*sqrt(1-0.5^2)*c/1e3;
    
    disp('Matrix Solved');
end
% Higher Order element
if abs(EleType - 2) < 1E-6
    disp('Building HD Influence Coefficient Matrix');
    disp('           ');
    %初始状态假设所有Joints上剪切应力方向为正
    [extraBC,CM] = BuildCoefMatix_H2(Mat,Fractures,AllEle_global(1:nAllEle_global,:),ConnList_global(1:nAllEle_global,:),1);
    %Unkown Displacements
    %Solve Linear System for Displacements
    disp('           ');
    disp('Matrix Solving...');
    disp('           ');
    BC = BC + extraBC;
    DD = CM\BC;
    disp('Matrix Solved');
end
for i = 1 : nAct
    DD_global(IndexInv(i)) = DD(i);
    DD_global(IndexInv(i) + MaxEle) = DD(i+nAct);
end
end