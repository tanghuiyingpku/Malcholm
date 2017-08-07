function [Vz,Vx] = CalcProppantV(nAct,Cp,Xf,dens_fl,propdens,visco_fl,w,Vel_sl)
global Fluid 
nprop = Fluid.nprop;
diam = zeros(nprop,1);
for ii = 1 : nprop
    diam(ii) = Fluid.prop{ii}.diam;
end
% kg/m^3 -> g/cm^3
propdens = propdens / 1e3;
% m -> cm
diam = diam * 100;
w = w * 100;

gee =980;
%pa.s --> p
visco_fl = visco_fl*10;
Vz = zeros(nAct*nprop,1);
Vx = zeros(nAct*nprop,5);
for i = 1 : nAct
    sumprop = sum(Cp(i,:));
    % Due to proppant concentration
    gc = 2.37*sumprop^2 - 3.08*sumprop + 1;
    for j = 1 : nprop
        denfl = dot(dens_fl(i,:),Xf(i,:));
        % g/cm^3
        denfl = denfl / 1000;
        % Stockes equilibrium velocity
        vstokes = (propdens(j) - denfl)*gee*diam(j)^2/18/visco_fl(i);
        % Due to turbulance
        vNe = 1.13*visco_fl(i)^0.57/denfl^0.29/(propdens(j) - denfl)^0.29/diam(j)^0.86;
        if diam(j) >= w(i)/2
            hw = 1;
        else
            wc = w(i);
            hw = 0.564*(diam(j)/wc)^2 - 1.563*diam(j)/wc + 1;
        end
        if Cp(i,j) > 1e-9
            %from cm/s to m/s
            Vz(i+(j-1)*nAct,1)  = vstokes*vNe*gc*hw/1e3;
            if diam(j) > w(i) || sumprop < 1e-6
                wc = diam(j)*1.6;
            else
                wc = 1/(1.411*(1/diam(j)^2 - 1/w(i)^2)*sumprop^0.8);
                wc = sqrt(wc);
            end
            k = 1+(diam(j)/wc) - 2.02*(diam(j)/wc)^2;
            % k*Vel_sl(i,:);
        end
        k=1;
        if Cp(i,j) > 1e-9
%             if diam(j)*2 > w(i)
%                 k=0;
%             else
%                 if diam(j)*4 > w(i)
%                     k = (w(i)-diam(j)*2)/2/diam(j);
%                 else
%                     k=1;
%                 end
%             end
        end
        Vx(i+(j-1)*nAct,:) =Vel_sl(i,:)*k;
    end
end