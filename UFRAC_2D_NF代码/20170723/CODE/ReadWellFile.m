function ReadWellFile(wellInput)
global nwell well Mat;
fid = fopen(wellInput,'r');
temp = fscanf(fid,'%s',1');
while ~feof(fid)
    if strcmp(temp, 'NofWells') == 1
        nwell = fscanf(fid,'%f',1);
        well = cell(nwell,1);
    end
    if strcmp(temp,'WellName') == 1
        for i = 1 : nwell
            well{i}.Name = fscanf(fid,'%s',1);
            well{i}.nSch = 0;
        end
    end
    if strcmp(temp, 'WellLocation') == 1
        for i = 1 : nwell
            well{i}.heel = fscanf(fid,'%f  %f',2);
            well{i}.toe = fscanf(fid,'%f  %f',2);
            well{i}.PerfP = -Mat.Sxx;
        end
    end
    if strcmp(temp, 'WellGridSize') == 1
        for i = 1 : nwell
            well{i}.nGrid = fscanf(fid,'%d',1);
        end
    end
    if strcmp(temp, 'Perfs') == 1
        for i = 1 : nwell
            temp = fscanf(fid,'%s',1);
            if strcmp(temp,well{i}.Name) == 1
                well{i}.nPerf = fscanf(fid,'%d',1);
                well{i}.Perf = zeros(well{i}.nPerf,2);
                well{i}.PerfType = zeros(well{i}.nPerf,1);
                for ii = 1 : well{i}.nPerf
                    well{i}.PerfType(ii) = fscanf(fid,'%d',1);
                end
                well{i}.PerfQ(1:well{i}.nPerf) = -999;
                for ii = 1 : well{i}.nPerf
                    well{i}.Perf(ii,1:2) = fscanf(fid,'%f %f',2);
                end
            end
        end
    end
    if strcmp(temp, 'Schedule') == 1 
        temp = fscanf(fid,'%s',1);
        while strcmp(temp,'END') < 0.1
            for i = 1 : nwell
                if strcmp(temp,well{i}.Name) == 1
                    well{i}.nSch = well{i}.nSch + 1;
                    well{i}.Sch(well{i}.nSch).t0 = fscanf(fid,'%f',1);
                    %Minutes to seconds
                    well{i}.Sch(well{i}.nSch).t0 = well{i}.Sch(well{i}.nSch).t0 * 60;
                    well{i}.Sch(well{i}.nSch).t1 = fscanf(fid,'%f',1);
                    well{i}.Sch(well{i}.nSch).t1 = well{i}.Sch(well{i}.nSch).t1 * 60;
                    well{i}.Sch(well{i}.nSch).wellType = fscanf(fid,'%s',1);
                    well{i}.Sch(well{i}.nSch).Contr = fscanf(fid,'%s',1);
                    well{i}.Sch(well{i}.nSch).ContrValue = fscanf(fid,'%f',1);
                    well{i}.Sch(well{i}.nSch).Prop = fscanf(fid,'%s',1);
                    well{i}.Sch(well{i}.nSch).PropFraction = fscanf(fid,'%f',1);
                    well{i}.Sch(well{i}.nSch).Fluid = fscanf(fid,'%s',1);
                    npf = fscanf(fid,'%d',1);
                    well{i}.Sch(well{i}.nSch).nPf = npf;
                    well{i}.Sch(well{i}.nSch).Pf = zeros(npf,1);
                    well{i}.Sch(well{i}.nSch).Pf_rate = zeros(npf,2);
                    well{i}.Sch(well{i}.nSch).Pf_Q = zeros(npf,2);
                    well{i}.Sch(well{i}.nSch).PfQsl = zeros(npf,2);
                    well{i}.Sch(well{i}.nSch).Pres =  zeros(npf,1);
                    for jj = 1 : npf
                        well{i}.Sch(well{i}.nSch).Pf(jj) = fscanf(fid,'%d',1);
                    end
                end
            end
            temp = fscanf(fid,'%s',1);
        end
    end
    temp = fscanf(fid,'%s',1');
end
WriteFracFile();

fclose(fid);
end