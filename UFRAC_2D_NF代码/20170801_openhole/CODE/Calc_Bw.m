function b = Calc_Bw(P,type)
%Punit:Mpa
%Water Compressibility
dP =5;
Pref = 1e5;% Unit:Pa
Bref = 1.00;
% From Mpa to Psi
%Bref = 1.00341;
%Incompressible
Cw = 0;%5*1e-8;
%Cw = 1;
Bw = Bref*exp((Pref - P)*Cw);
Bw2 = Bref*exp((Pref - (P+dP))*Cw);
if type == 1
    b = 1/Bw;
end
if type == 2
    b = (1/Bw2-1/Bw)/dP;
end
end