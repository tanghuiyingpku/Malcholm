function Padd = FindEquilP(Pinit,MassAtTip,xtip,Pbc,D,dx)
% nele = length(xtip);
% dmass = 1e6;
% Pt = Pinit;
% while abs(dmass) > 1e-5
%     Mass = 0;
%     for i = 1 : nele
%         pi = Pt/xtip(1)*xtip(i)+(1-xtip(i)/xtip(1))*Pbc;
%         b = Calc_Bw(pi,1);
%         Mass = Mass +  b*D(i)*dx(i);
%     end
%     dmass = (Mass-MassAtTip)/MassAtTip;
%     if dmass > 1e-6
%         Pt = Pt*(1-0.1);
%     else
%         Pt = Pt*(1+0.1);
%     end
%     Padd = Pt/xtip(1)*xtip+(1-xtip/xtip(1))*Pbc;
% end
Padd = (Pinit- Pbc)/xtip(1)*xtip+Pbc;
%Padd = Pt/xtip(1)*xtip+(1-xtip/xtip(1))*Pbc;
end
