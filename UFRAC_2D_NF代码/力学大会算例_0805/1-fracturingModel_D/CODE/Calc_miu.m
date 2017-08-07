function b = Calc_miu(P,type)
%Punit:Mpa
global Mat Fluid;

% viscosity Compressibility
dP = 5;
Pref = 1e5;%Bar surface pressure %psi to Mpa
Cmiu =0;% 1e-6;
miu_ref = Fluid.fluid{1}.visco;%Mat.visco;%100*1e-3;%;*100;% 100 cp
miu = miu_ref*exp((P - Pref)*Cmiu);
miu2 = miu_ref*exp((P+dP/1e6 - Pref)*Cmiu);
if type == 1
    b = 1/miu;
end
if type ==2
    b = (1/miu2 - 1/miu)/dP;
end
end