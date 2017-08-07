function SearchRadius()
global SearchRadius nFracture Fractures nAct  nActTip Tipcoordinate TipStatesInv;
global AllEle_global IndexInv Index ArrivalNF Mat XfF_global PresF_global ;
SearchRadius = 5;
for it = 1 : nActTip
    num = TipStatesInv(it);
    mx = Tipcoordinate(num,1);
    my = Tipcoordinate(num,2);
    Ele = AllEle_global(num,1:4);
    dir = [mx - AllEle_global(num,8),my - AllEle_global(num,9)];
    for i = 1 : nFracture
        % Calculate Distance
        if Fractures{i}.Type == 2
            Loc = Fractures{i}.Loc;
            [~, coord] = FindCrossPoint(Ele, Loc);
            isIn = isDotIn(coord,Loc);
            point = coord;
            if isIn > 0.1
                dis = CalculateDis(point,[mx,my]);
                dir1 = [point(1)-mx,point(2)-my];
                if dot(dir,dir1) > 1e-6
                    if dis < SearchRadius && ArrivalNF(i) < 0.1 
                        %Active This NF
                        startI = Fractures{i}.index(1);
                        ArrivalNF(i)=1;
                        for j = 1 : sum(Fractures{i}.eleNum)
                            nAct = nAct + 1;
                            IndexInv(nAct) = startI + j-1;
                            Index(startI+j-1) = nAct;
                            AllEle_global(startI + j-1,10) = 2;
                            PresF_global(startI + j-1) = Mat.Pp;
                            XfF_global(startI + j-1,:) = [1,0];
                        end
                    end
                end
            end
        end
    end
end
% uPDATE dd
% if hasChange > 0.1
%      CalcDD_nf();
% end
end