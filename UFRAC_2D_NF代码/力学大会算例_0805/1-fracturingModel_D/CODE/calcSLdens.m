function [den_sl,dens_dp] = calcSLdens(Pi,nFluid,fdens,fcmp,frefp,xfluid,cprop,propdens)
densf = zeros(nFluid,1);
den_fl = 0;
cprops = sum(cprop);
%fcmp = fcmp;%UNIT Mpa-1
for ii = 1 : nFluid
    densf(ii) = fdens(ii)*(1+fcmp(ii)*(Pi-frefp(ii)));
    den_fl = den_fl + (1-cprops)*xfluid(ii)*densf(ii);
end
den_sl = den_fl + dot(cprop,propdens);
dP = Pi*0.001;
Pi = Pi + dP;
den_fl = 0;
for ii = 1 : nFluid
    densf(ii) = fdens(ii)*(1+fcmp(ii)*(Pi-frefp(ii)));
    den_fl = den_fl + (1-cprops)*xfluid(ii)*densf(ii);
end
den_sl2 = den_fl + dot(cprop,propdens);
dens_dp = (den_sl2 - den_sl)/dP;
end