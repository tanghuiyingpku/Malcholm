function saveWellData(CurT)
global nwell well IndexInv stepL
global  MaxEle  AllEle_global
global    PresF_global nAct DD_global
global FILEPATH_main  Kindex RecordK Mat;
Path = ([FILEPATH_main,'\Data\']);
K_IC  =  RecordK(1:Kindex);
L = 0;
for i = 1 : nAct
    L = L + AllEle_global(IndexInv(i),7);
end
L = L/2 - stepL*2;  
fidL = fopen([Path,'HalfLength','.in'],'a+');
fprintf(fidL,'%f  %f\n', CurT,L);
fclose(fidL);
%mkdir([FILEPATH_main,'\Data']);

save KIC.txt -ascii K_IC;
copyfile('KIC.txt',[Path,'KIC.txt']);
for welli = 1 : nwell
    fidP = fopen([Path,'WellPresure',num2str(welli),'.in'],'a+');
    fidW = fopen([Path,'WellWidth',num2str(welli),'.in'],'a+');
    fidQ = fopen([Path,'WellQ',num2str(welli),'.in'],'a+');
    for i = 1 : well{welli}.nSch
        if CurT >= well{welli}.Sch(i).t0 && CurT <= well{welli}.Sch(i).t1
            Schindex = i;
        end
    end
    nPerf  = well{welli}.Sch(Schindex).nPf;
    Perf = well{welli}.Sch(Schindex).Pf;
    for i = 1 : nPerf
        numPf = Perf(i);
        numFrac = well{welli}.Perfindex(numPf,:);
        QQ = well{welli}.Sch(Schindex).Pf_Q(i,1); 
        Open = -DD_global(numFrac(1)+MaxEle);
        Open2 = -DD_global(numFrac(1)+MaxEle-1);
        Open = (3*Open - Open2)/2;
        Pres = PresF_global(numFrac(1));
        Pres2 = PresF_global(numFrac(1)-1);
        Pres = (3*Pres - Pres2)/2;
        fprintf(fidP,'%f  %d  %f\t\n',CurT, numPf,Pres);
        fprintf(fidW,'%f  %d  %f\t\n',CurT, numPf,Open);
        fprintf(fidQ,'%f  %d  %f\t\n',CurT, numPf,QQ);
    end
    fclose(fidP);
    fclose(fidW);
    fclose(fidQ);
end
end
