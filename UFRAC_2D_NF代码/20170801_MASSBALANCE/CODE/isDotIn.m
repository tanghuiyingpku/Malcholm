function [sta] = isDotIn(dot, L)
sta = 0;
V1 = dot - L(1:2);
V2 = dot - L(3:4);
% V3 = L(1:2)- L(3:4);
% D1 = norm(V1);
% D2 = norm(V2);
% D3 = norm(V3);
% f = abs(D1 + D2 - D3)/D3;
% if (f < eps)
%     sta = 1;
% end
cosdot = V1*V2';
if cosdot <= 1e-3
    sta = 1;
end
end
