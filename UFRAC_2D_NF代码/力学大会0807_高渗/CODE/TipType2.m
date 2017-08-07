function [type,vec] = TipType2(ne,AllEle,ConnList)
L = ConnList(ne,3);
R = ConnList(ne,4);
p0 = [AllEle(ne,8),AllEle(ne,9)];
p1 = [-999,-999];
if L > -0.1
    nbF = L;
else
    nbF = R;
end
if ConnList(ne,2) < 2.1
    p1 = [AllEle(nbF,8),AllEle(nbF,9)];
    if ConnList(nbF,1) ~= ConnList(ne,1)
        p1 = [AllEle(ne,8),AllEle(ne,9)];
    end
else
    M = ConnList(ne,5);
    if L > 0.1
        if ConnList(L,1) == ConnList(ne,1)
            p1 = [AllEle(L,8),AllEle(L,9)];
        else
            if R > 0.1
                if ConnList(R,1) == ConnList(ne,1)
                    p1 = [AllEle(R,8),AllEle(R,9)];
                else
                    if ConnList(M,1) == ConnList(ne,1)
                        p1 = [AllEle(M,8),AllEle(M,9)];
                    end
                end
            end
        end
    end
end
if p1(1) < - 800
    f = 1;
end
dis = p0-p1;
vec = dis/sqrt(dis(1)^2+dis(2)^2);
vec0 = [AllEle(ne,6) AllEle(ne,5)];
if dot(vec,vec0) < -1e-8
    %ÄæÏò
    type = 2;
else
    type = 1;
end

end