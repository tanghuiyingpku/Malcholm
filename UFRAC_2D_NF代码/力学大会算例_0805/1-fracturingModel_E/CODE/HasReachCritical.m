function [hasR,hasR1,nf1] = HasReachCritical(CritK)
global nActTip TipStates TipStatesInv;
hasR = zeros(50,1);
num = 0;
nf1 = 0;
% TipStates = 999 means the intersection location on HF
for i = 1 : nActTip
    nf = TipStatesInv(i);
    if TipStates(nf) > 998
        KC = CritK(i);
        num = num+ 1;
        if KC < 1 - 0.1
            hasR(num) = 0;
        else
            hasR(num) = 1 ;
            nf1 = nf;
        end
    end
end
if max(hasR) > 0.1
    hasR1 = 1;
else
    hasR1 = 0;
end
end