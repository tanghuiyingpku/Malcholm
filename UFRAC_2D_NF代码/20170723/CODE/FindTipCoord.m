function coord = FindTipCoord(ne)
global AllEle_global ConnList_global
nc = ConnList_global(ne,2);
point1 = AllEle_global(ne,1:2);
point2 = AllEle_global(ne,3:4);
in1 =0;
in2 = 0;
for i = 1 : nc
    nb = ConnList_global(ne,2+i);
    if nb > 0.1
        isIn = isDotIn(point1,AllEle_global(nb,1:4));
        isIn2 = isDotIn(point2,AllEle_global(nb,1:4));
        in1  = in1 + isIn;
        in2 = in2 + isIn2;
    end
end
if in1 < 0.1
    coord = point1;
else
    coord = point2;
end

end