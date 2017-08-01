function [w2,dw] = calcNfWidth(Sn,P)
global  InitialAperture
w0 = InitialAperture*5;
kal = 1e-8*5;
wr = InitialAperture*1;
% if P < P0
%     w2 = w0;
%     dw = 0;
% else
w2 = w0*exp(kal*(P-Sn))+wr;
dw = w0*kal*exp(kal*(P-Sn));
% end
% if P < -1e-6
%     w2 = InitialAperture;
%     dw = 0;
% end

if w2 > 0.1
    w2 = 0.1;
    dw = w0*kal*exp(kal*(log(w2/w0)/kal));
end

% global  InitialAperture
% w0 = InitialAperture*10;
% kal = 1e-9;
% % if P < P0
% %     w2 = w0;
% %     dw = 0;
% % else
% w2 = w0*exp(kal*(P-Sn));
% dw = w0*kal*exp(kal*(P-Sn));
% % end
% if P < -1e-6
%     w2 = InitialAperture;
%     dw = 0;
% end
% 
% if w2 > 0.1
%     w2 = 0.1;
%     dw = w0*kal*exp(kal*(log(w2/w0)/kal));
% end

end