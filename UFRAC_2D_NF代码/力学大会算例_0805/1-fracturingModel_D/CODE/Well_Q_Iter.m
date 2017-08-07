function dQ = Well_Q_Iter(P,welli,CurT)
% Using WU KAN' MODEL BUT WITH NewTonian Fluid
global well  Mat  Index
global Fluid  ActEle

for i = 1 : well{welli}.nSch
    if CurT >= well{welli}.Sch(i).t0 && CurT <= well{welli}.Sch(i).t1
        Schindex = i;
    end
end

%unknowns P0 Qf1,2,3...

%[P0 QF1... P0(Well2) QF1..]
% Pressure equation correspond to Q; Q equation correspond to P0

n = 2; % two sides
wellD = 0.1; % Well Diameter
nFluid = Fluid.nfluid;
nProp = Fluid.nprop;
proptype =  well{welli}.Sch(Schindex).Prop;
fluidtype = well{welli}.Sch(Schindex).Fluid;
nPerf  = well{welli}.Sch(Schindex).nPf;
Perf = well{welli}.Sch(Schindex).Pf;
fdens = zeros(nFluid,1);
frefp = zeros(nFluid,1);
fcmp = zeros(nFluid,1);
propdens = zeros(nProp,1);
fvisco = zeros(nFluid,1);
for i = 1 : nFluid
    fvisco(i) = Fluid.fluid{i}.visco;
    name = Fluid.fluid{i}.name;
    fdens(i) = Fluid.fluid{i}.dens;
    frefp(i) = Fluid.fluid{i}.refp;
    fcmp(i) = Fluid.fluid{i}.compr;
    if strcmp(name,fluidtype) == 1
        fluidindex = i;
    end
end
%
for i = 1 : nProp
    name = Fluid.prop{i}.name;
    propdens(i) = Fluid.prop{i}.dens;
    if strcmp(name,proptype) == 1
        propindex = i;
    end
end

xfluid = zeros(nFluid,1);
xfluid(fluidindex) = 1;
cprop = zeros(nProp,1);
cprop(propindex) =  well{welli}.Sch(Schindex).PropFraction;
isInj = 0;
if strcmp(well{welli}.Sch(Schindex).wellType,'INJ') == 1
    isInj = 1;
end
if strcmp(well{welli}.Sch(Schindex).wellType,'PROD') == 1
    isInj = -1;
end
Q0 = isInj * well{welli}.Sch(Schindex).ContrValue*0.159/60;
[dens_sl0,~] = calcSLdens(Mat.Pp*1e6,nFluid,fdens,fcmp,frefp,xfluid,cprop,propdens);
Dens = dens_sl0;
nVari = 0;
cmax = 0.7;
n = 1.3;
Vari = zeros(nPerf,1);
for i = 1 : nPerf
    numPf = Perf(i);
    numFrac = well{welli}.Perfindex(numPf,:);
    well{welli}.Sch(Schindex).PfQsl(i,1) = 0;
    well{welli}.Sch(Schindex).Pf_rate(i,1) =  0;
    well{welli}.Sch(Schindex).Pf_Q(i,1) =  0;
    well{welli}.Sch(Schindex).PfQsl(i,2) = 0;
    well{welli}.Sch(Schindex).Pf_rate(i,2) = 0;
    well{welli}.Sch(Schindex).Pf_Q(i,2) = 0;
    if ActEle(numFrac(1)) > 0.1
        nVari = nVari + 1;
        Vari(i) = i;
    end
    for mm = 1 : nProp
        well{welli}.Sch(Schindex).PfCp(i,mm) = cprop(mm);
    end
    for mm = 1 : nFluid
        well{welli}.Sch(Schindex).PfXf(i,mm) = xfluid(mm);
    end
end
cprops = sum(cprop);
visco_flv= dot(xfluid,fvisco);
visco_slv = visco_flv*(1-cprops/cmax)^-n;
nVari = nVari+ 1;

Jacob = zeros(nVari, nVari);
RHS = zeros(nVari,1);
x = zeros(nVari,1);
Qchara = 1;%Q0;
Pchara = 1;%abs(Mat.Sxx)*1e6;
Kd = 0.8;
Dperf = 20*1e-3*100/2.54;
Dens = Dens*0.0083454045;
n2 = 12;
% unit dens: lbm/gal dp:in RATE bpm = 0.159/60 m^3/s
Qunit = (1*60/0.159)^2;
K1 = 0.2369*Dens/n2^2/Dperf^4/Kd^2*Qchara^2/Pchara*Qunit;
% From psi to Pa
K1 = K1 / 145 * 1E6;
K1v = K1 + zeros(nVari-1,1);
K2 = 128*visco_slv/pi/wellD^4*Qchara/Pchara*1;

x(2:nVari) = Q0/(nVari-1)/Qchara;
x(1) = -Mat.Sxx*1.2*1e6/Pchara;
dis = 1e3;
tol = 1e-3;
nIter = 0 ;

while dis > tol
    nIter = nIter+1;
    RHS = RHS * 0;
    Jacob = Jacob *0;
    L = zeros(nVari-1,1);
    RHS(1) = sum(x(2:nVari)) - Q0;
    Jacob(1,2:nVari) = 1;
    for ip = 1 : nVari-1
        K1i = K1v(ip);
        RHS(ip+1) = x(1) -  P(Index(well{welli}.Perfindex(Vari(ip))))/Pchara - K1i*x(ip+1)^2 ;
        Jacob(ip+1,ip+1) = Jacob(ip+1,ip+1) - 2*K1i*x(ip+1);
        Jacob(ip+1,1) = 1;
        if ip == 1
            L(ip) = CalculateDis(well{welli}.heel', well{welli}.Perf(ip,:));
            dPfric = K2*L(ip)*Q0/Qchara;
        else
            L(ip)  = CalculateDis(well{welli}.Perf(ip,:), well{welli}.Perf(ip-1,:));
            dPfric = K2*L(1)*Q0/Qchara;
            for jj = 1: ip-1
                dPfric = dPfric + L(jj+1)*(Q0/Qchara - sum(x(2:jj+1)))* K2;
                Jacob(ip+1,jj+1) = Jacob(ip+1,jj+1) + sum(L(jj+1:ip)) * K2;
            end
        end
        RHS(ip+1) = RHS(ip+1) - dPfric ;
    end
    dx = -Jacob\RHS;
    x = x + dx;
    dis = norm(RHS);
end
addQ = 0;
num = 0;
for i = 1 : nVari
    if x(i) < 1e-10
        num = num +1;
        addQ = addQ-x(i);
        x(i) = 1e-5;
    end
end
nEffec = nVari-1-num;
aveQ = addQ/nEffec;
for i = 2 : nVari
    if x(i) > 1e-4
        x(i) = x(i)-aveQ;
    end
end
x(1) = x(1) * Pchara;
x(2:nVari) = x(2:nVari) * Qchara ;
Q2 = x(2:nVari);
dQ = 0;
for i = 1 : nVari-1
    well{welli}.Sch(Schindex).PfQsl(Vari(i),1) = dens_sl0*Q2(i)/2;
    well{welli}.Sch(Schindex).PfQsl(Vari(i),2) = dens_sl0*Q2(i)/2;
    dQ = abs(well{welli}.Sch(Schindex).Pf_Q(Vari(i),1) - Q2(i)/2)+dQ;
    well{welli}.Sch(Schindex).Pf_Q(Vari(i),1) = Q2(i)/2;
    well{welli}.Sch(Schindex).Pf_Q(Vari(i),2) = Q2(i)/2;
end
%well{welli}.P = x(1);
%        well{welli}.P = x(1);
end