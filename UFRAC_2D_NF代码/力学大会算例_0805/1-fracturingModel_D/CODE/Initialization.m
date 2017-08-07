function [DD,Pres] = Initialization(nFracture,Fractures,nAllEle,AllEle,ConnList,Mat,BC,EleType,PP)
global Index IndexInv nAct
global MaxEle;
IndexInv = zeros(MaxEle,1);
Index = zeros(MaxEle,1);
%BC 部分不用修改 已经生成好
%EleType 表示使用高阶还是非高阶进行计算
disp('           ');
disp('***************************************************************');
disp('Initialization Begin: Calculate the initial stress distribution');
disp('           ');
disp('Height Correction(Olson 2004) is used');
disp('           ');
cnt = 1;
for ii = 1 : nAllEle
    if (abs(AllEle(ii,10) - 3) < 1e-6) ||  abs(AllEle(ii,10) - 1) < 1e-6
        IndexInv(cnt) = ii;
        Index(ii) = cnt;
        cnt = cnt + 1;
    end
end
nAct = cnt-1;
if abs(EleType - 1) < 1E-6
    disp('Building CCD Influence Coefficient Matrix');
    disp('           ');
    [extraBC,CM]  = BuildCoefMatix_Constant(Mat,nFracture,Fractures,nAllEle,AllEle,ConnList,1);
    %Unkown Displacements
    %Solve Linear System for Displacements
    disp('           ');
    disp('Matrix Solving...');
    disp('           ');
    
    BCi = zeros(cnt-1,1);
    Pres = zeros(cnt-1,1);
    for ii = 1 : cnt-1
        Pres(ii) = PP(IndexInv(ii));
        BCi(ii) = BC(IndexInv(ii));
        BCi(ii+cnt-1) = BC(IndexInv(ii)+nAllEle);
    end
    BCi = BCi*1e6 + extraBC;
    DD = CM\BCi;
    disp('Matrix Solved');
end
% Higher Order element
if abs(EleType - 2) < 1E-6
    disp('Building HD Influence Coefficient Matrix');
    disp('           ');
    %初始状态假设所有Joints上剪切应力方向为正
    [extraBC,CM] = BuildCoefMatix_H2(Mat,nFracture,Fractures,nAllEle,AllEle,ConnList,1);
    %Unkown Displacements
    %Solve Linear System for Displacements
    disp('           ');
    disp('Matrix Solving...');
    disp('           ');
    BCi = zeros(cnt-1,1);
    Pres = zeros(cnt-1,1);
    for ii = 1 : cnt-1
        Pres(ii) = PP(IndexInv(ii));
        BCi(ii) = BC(IndexInv(ii));
        BCi(ii+cnt-1) = BC(IndexInv(ii)+nAllEle);
    end
    BCi = BCi*1e6 + extraBC;
    DD = CM\BCi;
    disp('Matrix Solved');
    
end
end