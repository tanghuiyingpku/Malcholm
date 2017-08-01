function WriteFracFile
global nwell well stepL
global Mat
global  FILEPATH_main
global nInitial;
nInitial = 3;
% nInitial
Ks = 10;
Kn = 10;
Fang = 20;
Coh = 1;

maxNfrac = 30;
frac = zeros(maxNfrac,4);
fracnf = zeros(maxNfrac,4);
fracFile =[FILEPATH_main, 'Input_frac.in'];
nffracFile =[FILEPATH_main, 'Input_nf.in'];
kinkfracFile =[FILEPATH_main, 'Input_kink.in'];
kinkfrac = load(kinkfracFile);
[n_kink,~] = size(kinkfrac);
nffrac = load(nffracFile);
nf1 = 0;
[n_nf,~] = size(nffrac);
L = stepL*nInitial;

% Natural Frac Part
for i = 1 : n_nf
    Li  = CalculateDis(nffrac(i,1:2),nffrac(i,3:4));
    ni = 2;%floor(Li/L);
    xx= linspace(nffrac(i,1),nffrac(i,3),ni);
    yy= linspace(nffrac(i,2),nffrac(i,4),ni);
    fracnf(nf1+1:nf1+ni-1,1) = xx(1:ni-1);
    fracnf(nf1+1:nf1+ni-1,2) = yy(1:ni-1);
    fracnf(nf1+1:nf1+ni-1,3) = xx(2:ni);
    fracnf(nf1+1:nf1+ni-1,4) = yy(2:ni);
    nf1 = nf1 + ni-1;
end
fid = fopen(fracFile,'w');
dir = Mat.SminDir+90;
%Matdir 0~90 degree
if dir > 90
    dir = dir - 180;
end
vec = [cosd(dir),sind(dir)];
nfrac = 0;
% Frac = nPerf
for i = 1 : nwell
    nperf = well{i}.nPerf;
    for j = 1 : nperf
        if well{i}.PerfType(j) > 0.1;
            center = well{i}.Perf(j,1:2);
            frac(nfrac+1,:) = [center-vec*L,center+vec*L];
            nfrac = nfrac + 1;
        end
    end
end
nfracII = nfrac;
for i = 1 : n_kink
    frac(nfrac+1,:) = kinkfrac(i,:);
    nfrac = nfrac + 1;
end
fprintf(fid,'-- number of fracture\n');
fprintf(fid,'%d\n',nfrac+nf1);
fprintf(fid,'-- fracture type 1:hydraulic fracture 2: joint\n');
for i = 1 : nfrac
    
    if i > nfracII
        fprintf(fid,'%d\n',7);
    else
        fprintf(fid,'%d\n',1);
    end
end
for i = 1 : nf1
    fprintf(fid,'%d\n',2);
end
fprintf(fid,'--  Frac element initial length of boundary elements\n');
for i = 1 : nfrac
    fprintf(fid,'%f\n',stepL);
end
for i = 1 : nf1
    fprintf(fid,'%f\n',stepL);
end



fprintf(fid,'--  locations of Frac start and end\n');
for i = 1 : nfrac
    fprintf(fid,'%f  %f  %f  %f\n',frac(i,:));
end
for i = 1 : nf1
    fprintf(fid,'%f  %f  %f  %f\n',fracnf(i,:));
end
fprintf(fid,'--   Boundary Value(BC)\n');
for i = 1 : nfrac
    fprintf(fid,'%f\n%f\n',0,0);
end
for i = 1 : nf1
    fprintf(fid,'%f\n%f\n%f\n%f\n',Ks,Kn,Fang,Coh);
end
fclose(fid);
clear frac fracnf;

end