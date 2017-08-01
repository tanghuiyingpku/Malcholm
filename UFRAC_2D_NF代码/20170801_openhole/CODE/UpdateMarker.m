function UpdateMarker(nEle,Conn)
global CM_marker;
global GrowthN;
find = 0;
for i = 1 : nEle
    if abs(CM_marker(i) - 3)<1e-3
        %最后一个编号为3的是Tip
        if find < 0.1
            start = i;
        end
        find = find + 1;
    end
    if abs(find-GrowthN) < 1e-3
        %get to the tip element
        if CM_marker(Conn(start,3)) < 2.9
            CM_marker(Conn(start,3)) = -1;
            nb = Conn(start,3);
            CM_marker(Conn(nb,3)) = -1;
            CM_marker(Conn(nb,4)) = -1;
            L = Conn(nb,3);
            R = Conn(nb,4);
            if CM_marker(Conn(L,3)) > -1e-3
                CM_marker(Conn(L,3)) = -1;
            end
            if CM_marker(Conn(L,4)) > -1e-3
                CM_marker(Conn(L,4)) = -1;
            end
            if CM_marker(Conn(R,3)) > -1e-3
                CM_marker(Conn(R,3)) = -1;
            end
            if CM_marker(Conn(R,4)) > -1e-3
                CM_marker(Conn(R,4)) = -1;
            end
        else
            CM_marker(Conn(start,4)) = -1;
            nb = Conn(start,4);
            CM_marker(Conn(nb,3)) = -1;
            CM_marker(Conn(nb,4)) = -1;
            L = Conn(nb,3);
            R = Conn(nb,4);
            if CM_marker(Conn(L,3)) > -1e-3
                CM_marker(Conn(L,3)) = -1;
            end
            if CM_marker(Conn(L,4)) > -1e-3
                CM_marker(Conn(L,4)) = -1;
            end
            if CM_marker(Conn(R,3)) > -1e-3
                CM_marker(Conn(R,3)) = -1;
            end
            if CM_marker(Conn(R,4)) > -1e-3
                CM_marker(Conn(R,4)) = -1;
            end
        end
        CM_marker(i) = 2;
        CM_marker(i-1) = 1;
        left = Conn(i-1,3);
        right = Conn(i-1,4);
        if  left ~= i
            CM_marker(left) = 0;
        else
            CM_marker(right) = 0;
        end
        find = 0;
    end
end
end