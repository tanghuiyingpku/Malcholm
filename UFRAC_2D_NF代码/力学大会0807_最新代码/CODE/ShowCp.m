function ShowCp(Cp)
dx = 1;
[nEle,nProp ] = size(Cp);
figure(11);
hi = zeros(nEle,1);
for i = 1 : nProp
    C = 90+(i-1)*50;
    for j = 1 : nEle
        x = [dx*(j-1),dx*j,dx*j,dx*(j-1)];
        y = [hi(j) , hi(j) , hi(j)+Cp(j,i), hi(j)+Cp(j,i)];
        patch(x,y,C);
        hold on
        hi(j) = hi(j) + Cp(j,i);
    end
end
ylim([0,1]);
ylabel('Proppant Distribution','Fontsize',12);
title('Proppant Distribution','Fontsize',13);
hold off;
end