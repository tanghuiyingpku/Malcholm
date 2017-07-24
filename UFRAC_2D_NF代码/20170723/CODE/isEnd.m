function isend = isEnd()
global nAct AllEle_global IndexInv
isend = 0;
for i = 1 : nAct
    if abs(AllEle_global(IndexInv(i),9)) > 200 || (AllEle_global(IndexInv(i),9)) < 0
        % if AllEle_global(IndexInv(i),9)< 0
        isend = 1;
    end
end