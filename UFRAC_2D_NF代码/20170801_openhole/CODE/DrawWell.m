function DrawWell(P,Xf,Cp)
global PicScale Fluid well;
n =  well{1}.nGrid;
figure(21);
plot(P,'.');
title('Pressure in Well','Fontsize',14);
figure(22);
color = ['b','m','c','r','g'];
for i = 1 : Fluid.nfluid;
    num = (i-1)*n;
    plot(Xf(num+1:num+n),[color(i),'.']);
    hold on;
end
for i = 1: Fluid.nprop
    num = (i-1)*n;
    plot(Cp(num+1:num+n),[color(i),'--']);
    hold on;
end
title('Components Distribution in Well','Fontsize',14);
ylim([-0.5 1.5]);
end