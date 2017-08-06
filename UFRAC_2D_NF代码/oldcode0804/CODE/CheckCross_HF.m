function [iHF,jHF,Point,isMerge,mindis] = CheckCross_HF(num,Ele,nAllEle,AllEle,Conn,L,nAct0)
global nAct IndexInv GrowthN Tipcoordinate;
[sizeC,~] = size(Conn);
eps = L;
mindis = 1e6;
iHF = -1;
jHF = -1;
HFnumEle0 = -1;
Point = [0 0];
isMerge = 0;
mx = Tipcoordinate(num,1);
my = Tipcoordinate(num,2);
dir = [mx - AllEle(num,8),my - AllEle(num,9)];
cp = zeros(2,1);
hfIndex_current = Conn(num,1);
for i = 1 : nAct0
    ic = IndexInv(i);
    if ic > sizeC
        f = 1;
    end
    if Conn(ic,1) ~= hfIndex_current && AllEle(ic,10) < 1.1
        Loc = AllEle(ic,1:4);
        [~, coord] = FindCrossPoint(Ele, Loc);
        isIn = isDotIn(coord,Loc);
        point = coord;
        if isIn > 0.1
           dis = CalculateDis(point,[mx,my]);
            dir1 = [point(1)-mx,point(2)-my];
            if dot(dir,dir1) > 1e-6
                if dis < mindis
                    mindis = dis;
                    HFnumEle0 = i;
                    HFnum = Conn(ic,1);
                    cp = point;
                end
            end
        end
    end
end 
if HFnumEle0 < 0.1
    return;
end
mindis2 = 1e6;
HFnumEle = -1;
for i = 1 : nAct-2
    ic = IndexInv(i);
    if Conn(ic,1) == HFnum
        dis = CalculateDis(cp,AllEle(ic,8:9));
        if dis < mindis2
            mindis2 = dis;
            HFnumEle = ic;
        end
    end
end
if HFnumEle < -0.1
    return;
end

if mindis < L
    %
    isMerge = 1;
    left = Conn(HFnumEle,3);
    right = Conn(HFnumEle,4);
    [~, coord] = FindCrossPoint(Ele, AllEle(HFnumEle,1:4));
    isIn = isDotIn(coord,AllEle(HFnumEle,1:4));
    if isIn > 0.1
        Point = coord;
        dis1 = CalculateDis(coord,AllEle(HFnumEle,1:2));
        dis2 = CalculateDis(coord,AllEle(HFnumEle,3:4));
        [~,b] = min([dis1,dis2]);
        if min(dis1,dis2) < eps
            %相交在两个element的交点上
            if b > 1.9
                jHF = right;
            else
                jHF = left;
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
            HFnumEle = left;
            jHF = HFnumEle;
        end
        if isIn > 0.1
            Point = coord;
            dis1 = CalculateDis(coord,AllEle(left,1:2));
            dis2 = CalculateDis(coord,AllEle(left,3:4));
            [~,b] = min([dis1,dis2]);
            if min(dis1,dis2) < eps
                %相交在两个element的交点上
                if b > 1.9
                    jHF = Conn(left,4);
                else
                    jHF = Conn(left,3);
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
            HFnumEle = right;
            jHF = HFnumEle;
        end
        if isIn > 0.1
            Point = coord;
            dis1 = CalculateDis(coord,AllEle(right,1:2));
            dis2 = CalculateDis(coord,AllEle(right,3:4));
            [~,b] = min([dis1,dis2]);
            if min(dis1,dis2) < eps
                %相交在两个element的交点上
                if b > 1.9
                    jHF = Conn(right,4);
                else
                    jHF = Conn(right,3);
                end
            end
        end
    end
    if isIn < 0.1
        isMerge = 0;
        return;
    end
end
iHF = HFnumEle;
end
