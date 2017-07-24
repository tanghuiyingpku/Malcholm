function ShowBanking(hbank)
global Mat;
dx = 1;
[nEle,nProp ] = size(hbank);
figure(10);
hi = zeros(nEle,1);
for i = 1 : nProp
    C = 90+(i-1)*50;
    for j = 1 : nEle
        x = [dx*(j-1),dx*j,dx*j,dx*(j-1)];
        y = [hi(j) , hi(j) , hi(j)+hbank(j,i), hi(j)+hbank(j,i)];
        patch(x,y,C);
        hold on
        hi(j) = hi(j) + hbank(j,i);
    end
end
ylim([0, max(hi)*3+1e-3]);
ylabel('Banking Height','Fontsize',12);
title('Banking Height','Fontsize',13);
hold off;

end