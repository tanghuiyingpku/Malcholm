function CompareAna()
global FILEPATH_main;
Path = ([FILEPATH_main,'\Data\']);
data = load( [Path,'HalfLength.in']);
t = data(:,1);
t = t - 0.2;
L = data(:,2);
[L2,W2,P2] = AnalyticalKGD(t);
data  = load([Path,'WellPresure1.in']);
pw = data(:,3);
data = load([Path,'WellWidth1.in']);
w = data(:,3);
figure(22)
plot(t/60,L,'b.','Linewidth',2);
hold on
plot(t/60,L2,'r--','Linewidth',2);
legend('Numerical','Analytical');
title('Half L');
hold off;

figure(23)
plot(t/60,w*100,'b.','Linewidth',2);
hold on
plot(t/60,W2*100,'r--','Linewidth',2);
legend('Numerical','Analytical');
title('Width');
hold off;

figure(24)
plot(t/60,pw,'b.','Linewidth',2);
hold on
plot(t/60,P2,'r--','Linewidth',2);
legend('Numerical','Analytical');
title('Pressure');
hold off;
end