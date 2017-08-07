function  FluidPressure = FluidEquation(Mat,P_pre,DD,dt,EleL,EleConn,ElementType)
% Unit system: Metric
% This part we use the simplest method
% FVM 
% and pure water fluid
% Assumptions: Isothermal T = 60 degree
% Pore pressure
h = 10;% 10m
Qinj = 3;% 5m^3/s
Pp = Mat.Pp*1e6;
P_pre = P_pre*1e6;

nEle = length(ElementType);
%viscosity of fluid
miu = 46.88*1e-5;
%Fluid Bulk modules
Ew = 1e3;

Pnew = P_pre;
Pold =  P_pre+Pp;

%Normal Opening
e = DD(nEle+1:2*nEle);
%Volumn
V = e.*EleL*h;
krock = 1e-15*1e-3;

% get dt time step
for i = 1 : nEle
    dti = 12*miu*EleL(i)^2/Ew/e(i)^2;
    if dti < dt
        dt = dti;
    end
end
%d is effective leak distance
d = 100;
dif = norm(Pnew-Pold);
eps = dif*1e-4;
MaxIter = 200;
iIter = 0;

while dif > eps
    iIter = iIter + 1;
    %Inflow
    Q = zeros(nEle,1);
    Qleak = zeros(nEle,1);
    P = Pnew;
    for i = 1: nEle
        %Flow
       con = EleConn(i,:);
       ki = e(i)^3/12/miu;
       Li = EleL(i);
       for j = 1 : length(con)
           if con(j) > -0.1
               kj =  e(con(j))^3/12/miu;
               Lj = EleL(con(j));
               kij = ki*kj*(Li+Lj)/(ki*Li+kj*Lj);
               Q(i) = Q(i) + kij*(P(con(j)) -P(i))*2/(Li+Lj);
           end
       end
       Q(i) = Q(i) * h;
       %Leak
       Qleak(i) = krock*Li/miu*(P(i)-Pp)/d*h;
       %d is effective leak distance
       if i == Mat.Perf
           Pnew(i) = P_pre(i) + Ew*Q(i)*dt/V(i) - Ew*Qleak(i)*dt/V(i)+ Ew*Qinj*dt/V(i);
       else
           Pnew(i) = P_pre(i) + Ew*Q(i)*dt/V(i) - Ew*Qleak(i)*dt/V(i);
       end
    end
       Pold = P;
       dif = norm(Pnew-Pold);
       if(iIter > MaxIter)
           error('Can not convergent within given iteration number');
       end
end
disp('Fluid and solid Coupling solved...');
fprintf('Iteration = %d\n',iIter);
FluidPressure = Pnew;
end