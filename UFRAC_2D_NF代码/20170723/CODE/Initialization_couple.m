function Initialization_couple()
global Index IndexInv nAct ActEle
global MaxEle;
IndexInv = zeros(MaxEle,1);
Index = zeros(MaxEle,1);
ActEle = zeros(MaxEle,1);
%BC 部分不用修改 已经生成好
%EleType 表示使用高阶还是非高阶进行计算
disp('           ');
disp('***************************************************************');
disp('Initialization Begin: Calculate the initial stress distribution');
disp('           ');
disp('Height Correction(Olson 2004) is used');
disp('           ');

nAct = 0;

end