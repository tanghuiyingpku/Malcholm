function Initialization_couple()
global Index IndexInv nAct ActEle
global MaxEle;
IndexInv = zeros(MaxEle,1);
Index = zeros(MaxEle,1);
ActEle = zeros(MaxEle,1);
%BC ���ֲ����޸� �Ѿ����ɺ�
%EleType ��ʾʹ�ø߽׻��ǷǸ߽׽��м���
disp('           ');
disp('***************************************************************');
disp('Initialization Begin: Calculate the initial stress distribution');
disp('           ');
disp('Height Correction(Olson 2004) is used');
disp('           ');

nAct = 0;

end