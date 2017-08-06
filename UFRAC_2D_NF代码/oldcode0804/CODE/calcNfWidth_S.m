function [w2,dw] = calcNfWidth_S(Ds,Sn,P,sigmaN)
global  InitialAperture Mat
%
Ds=0;
S0 = 20;
SN0 =0;%Mat.Pp*1e6 + sigmaN;
phi = 10/180*pi;
%
if abs(Ds)>1e-6
    f = 1;
end
w0 = InitialAperture*5;
kal = 1e-7;
wr = InitialAperture*1;
w2 = w0*exp(kal*(P-Sn-SN0))+ wr + abs(Ds)*tan(phi/(1+9*(Sn-P)/1e6/S0));
P2 = P + 10;
w22 = w0*exp(kal*(P2-Sn-SN0))+ wr + abs(Ds)*tan(phi/(1+9*(Sn-P2)/1e6/S0));
dw = (w22-w2)/10;
%

if w2 > 0.1
    w2 = 0.1;
    dw = w0*kal*exp(kal*(log(w2/w0)/kal));
end

% 
% S0 = 20;
% w0 = InitialAperture;
% wr = InitialAperture*0.5;
% phi = 0/180*pi;
% 
% if abs(Ds) > 1e-2
%     Ds = 1e-2;
% end
% 
% 
% w2 = w0/(1+9*(Sn-P)/1e6/S0) + abs(Ds)*tan(phi/(1+9*(Sn-P)/1e6/S0)) + wr;
% dP = P*1e-5;
% if abs(P) < 1e-5
%     dP = 1;
% end
% P2 = P + dP;
% w22 = w0/(1+9*(Sn-P2)/1e6/S0) + abs(Ds)*tan(phi/(1+9*(Sn-P2)/1e6/S0)) + wr;
% dw = (w22-w2)/dP;
% if abs(P) < 1e-5
%     dw=0;
% end

end