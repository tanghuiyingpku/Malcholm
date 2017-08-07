function     DrawJacob(Jac)
figure(10);
[nx,ny] = size(Jac);
% for i = 1 : nx
%     for j = 1 : ny
%         if abs(Jac(i,j))>1e3
%             Jac(i,j) = 1e-4;
%         end
%     end
% end
imagesc(Jac);
title('Jacob Structure','Fontsize',14);
end
