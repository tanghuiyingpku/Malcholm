function saveData()
global  AllEle_global  PicScale;
global  MaxEle nAllEle_global
global    PresF_global e_global DD_global
global FILEPATH_main nAct IndexInv;
mkdir([FILEPATH_main,'\Data']);
Path = ([FILEPATH_main,'\Data\']);
Element = zeros(nAllEle_global,5);
OpenN = zeros(nAllEle_global,1);
OpenS = zeros(nAllEle_global,1);
Pres = zeros(nAllEle_global,1);
HydrauE = zeros(nAllEle_global,1);
OpenA = zeros(nAct,1);
HydrauEA = zeros(nAct,1);
SA = zeros(nAct,1);
PA = zeros(nAct,1);
XA = zeros(nAct,1);
YA = zeros(nAct,1);

for i = 1 : nAct
    SA(i) = DD_global(IndexInv(i));
    PA(i) = PresF_global(IndexInv(i));
    XA(i) = AllEle_global(IndexInv(i),8);
    YA(i) = AllEle_global(IndexInv(i),9);
    OpenA(i) = -DD_global(IndexInv(i)+MaxEle);
    HydrauEA(i) = e_global(IndexInv(i));
end
% Marker HF 1 close NF -1 open NF 2
for i = 1 : nAllEle_global
    ele = AllEle_global(i,1:4);
    if abs(AllEle_global(i,10) - 1) < 0.1
        frac = [ele,1];
    end
     if abs(AllEle_global(i,10) - 2) < 0.1
        frac = [ele,-1];
     end
     if abs(AllEle_global(i,10) - 3) < 0.1
        frac = [ele,2];
     end
    ds = DD_global(i);
    dn = -DD_global(i+MaxEle);
    pi = PresF_global(i);
    e = e_global(i);
    Element(i,:) = frac;
    OpenN(i) = dn;
    OpenS(i) = ds;
    Pres(i) = pi;
    HydrauE(i) = e;
end

save Element.txt -ascii Element;    
save OpenN.txt -ascii OpenN;    
save OpenS.txt -ascii OpenS;    
save Pres.txt -ascii Pres;    
save HydrauE.txt -ascii HydrauE;   
save HydrauEA.txt -ascii HydrauEA;    
save OpenA.txt -ascii OpenA;    
save SA.txt -ascii SA;    
save PA.txt -ascii PA;    
save XA.txt -ascii XA;    
save YA.txt -ascii YA;    

save PicScale.txt -ascii PicScale;    
copyfile('Element.txt',[Path,'Element.txt']);
copyfile('OpenN.txt',[Path,'OpenN.txt']);
copyfile('OpenS.txt',[Path,'OpenS.txt']);
copyfile('OpenA.txt',[Path,'OpenA.txt']);
copyfile('SA.txt',[Path,'SA.txt']);
copyfile('XA.txt',[Path,'XA.txt']);
copyfile('YA.txt',[Path,'YA.txt']);
copyfile('PA.txt',[Path,'PA.txt']);

copyfile('Pres.txt',[Path,'Pres.txt']);
copyfile('HydrauE.txt',[Path,'HydrauE.txt']);
copyfile('HydrauEA.txt',[Path,'HydrauEA.txt']);
copyfile('PicScale.txt',[Path,'PicScale.txt']);

end