function [iNF,jNF,Point,isMerge,mindis] = CheckCross(num,Ele,nAllEle,AllEle,Conn,L)
global Fractures;
global nFracture Tipcoordinate;
NFnum = -1;
eps = L*0.5;
mindis = 1e6;
iNF = -1;
jNF = iNF;
Point = [0 0];
isMerge = 0;
mx = Tipcoordinate(num,1);
my = Tipcoordinate(num,2);
dir = [mx - AllEle(num,8),my - AllEle(num,9)];
cp = zeros(2,1);
for i = 1 : nFracture
    %如果是天然裂缝
    if Fractures{i}.Type == 2
        Loc = Fractures{i}.Loc;
        [~, coord] = FindCrossPoint(Ele, Loc);
        isIn = isDotIn(coord,Loc);
        point = coord;
        if isIn > 0.1
            dis = CalculateDis(point,[mx,my]);
            dir1 = [point(1)-mx,point(2)-my];
            if dot(dir,dir1) > 1e-6
                if dis < mindis
                    mindis = dis;
                    NFnum = i;
                    cp = point;
                end
            end
        end
    end
end
if NFnum < 0.1
    return;
end
mindis2 = 1e6;
for i = 1 : nAllEle
    if Conn(i,1) == NFnum
        dis = CalculateDis(cp,AllEle(i,8:9));
        if dis < mindis2
            mindis2 = dis;
            iNF = i;
        end
    end
end
if iNF < -0.1
    return;
end
%找到NF之后计算到该NF的端点连线距离
jNF = iNF;
if mindis < L
    %
    isMerge = 1;
    left = Conn(iNF,3);
    right = Conn(iNF,4);
    [~, coord] = FindCrossPoint(Ele, AllEle(iNF,1:4));
    isIn = isDotIn(coord,AllEle(iNF,1:4));
    if isIn > 0.1
        Point = coord;
        dis1 = CalculateDis(coord,AllEle(iNF,1:2));
        dis2 = CalculateDis(coord,AllEle(iNF,3:4));
        [~,b] = min([dis1,dis2]);
        if min(dis1,dis2) < eps
            %相交在两个element的交点上
            if b > 1.9
                jNF = right;
            else
                jNF = left;
            end
        end
    end
    
    if isIn < 0.1
        if left < -0.1
            coord = [-999,-999];
            isIn = 0;
        else
            [~, coord] = FindCrossPoint(Ele, AllEle(left,1:4));
            isIn = isDotIn(coord,AllEle(left,1:4));
            iNF = left;
            jNF = iNF;
        end
        if isIn > 0.1
            Point = coord;
            dis1 = CalculateDis(coord,AllEle(left,1:2));
            dis2 = CalculateDis(coord,AllEle(left,3:4));
            [~,b] = min([dis1,dis2]);
            if min(dis1,dis2) < eps
                %相交在两个element的交点上
                if b > 1.9
                    jNF = Conn(left,4);
                else
                    jNF = Conn(left,3);
                end
            end
        end
    end
    if isIn < 0.1
        if right < -0.1
            coord = [-999,-999];
            isIn = 0;
        else
            [~, coord] = FindCrossPoint(Ele, AllEle(right,1:4));
            isIn = isDotIn(coord,AllEle(right,1:4));
            iNF = right;
            jNF = iNF;
        end
        if isIn > 0.1
            Point = coord;
            dis1 = CalculateDis(coord,AllEle(right,1:2));
            dis2 = CalculateDis(coord,AllEle(right,3:4));
            [~,b] = min([dis1,dis2]);
            if min(dis1,dis2) < eps
                %相交在两个element的交点上
                if b > 1.9
                    jNF = Conn(right,4);
                else
                    jNF = Conn(right,3);
                end
            end
        end
    end
    if isIn < 0.1
        isMerge = 0;
        return;
    end
end

end
