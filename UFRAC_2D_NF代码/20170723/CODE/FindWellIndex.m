function FindWellIndex()
global nwell well nAllEle_global AllEle_global ConnList_global
target = 1;
for i = 1 : nwell
    nPerf = well{i}.nPerf;
    nG = well{i}.nGrid;
    well{i}.Gridx = linspace(well{i}.heel(1),well{i}.toe(1),nG+1);
    well{i}.Gridy = linspace(well{i}.heel(2),well{i}.toe(2),nG+1);
    well{i}.dx =  CalculateDis([well{i}.Gridx(1),well{i}.Gridy(1)],[well{i}.Gridx(2),well{i}.Gridy(2)]);
    for ii = 1 : nPerf
        loc = well{i}.Perf(ii,:);
        for mm = 1 : nG-1
            if loc(1) >= well{i}.Gridx(mm) && loc(1) <= well{i}.Gridx(mm+1)
                well{i}.PerfWellindex(ii) = mm;
                break;
            end
        end
        dismin = 1e6;
        for j = 1 : nAllEle_global
            if well{i}.PerfType(ii) > 0.1
                if AllEle_global(j,10) > 1.1
                    continue;
                end
            end
            xm = AllEle_global(j,8);
            ym = AllEle_global(j,9);
            dis =  CalculateDis([xm,ym],loc);
            if dis < dismin
                dismin = dis;
                target = j;
            end
        end
        tti = ConnList_global(target,3);
        xm = AllEle_global(tti,8);
        ym = AllEle_global(tti,9);
        dis =  CalculateDis([xm,ym],loc);
        if dis > dismin*1.1
            tti =  ConnList_global(target,4);
        end
        well{i}.Perfindex(ii,:) = sort([target;tti]);
    end
end
end