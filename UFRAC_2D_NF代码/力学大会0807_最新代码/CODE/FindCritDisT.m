function [Lleft2,Lright2,L2,hasFind] = FindCritDisT(Lleft,Lright,L,sigmaN,To,tol)
Lleft0 = 1e-5;
Lright0 = 1e-1;
Lright2 = Lright;
Lleft2 = Lleft;
L2 = L;
hasFind = 0;
coeff = -sigmaN+To;
if abs(coeff) > To*tol
    hasFind = 0;
    if coeff >= To*tol % KIC is large
        if Lleft <= Lleft0 ||  Lright >= Lright0
            Lright2 = L;
            L2 = L/1.5;
        end
        if Lright < Lright0 && Lleft > Lleft0
            L2 = 0.5*(Lleft+ Lright);
        end
        if Lright < Lright0 && L < Lright
            Lright2 = L;
        end
    end
    if coeff <= -To*tol
        if Lleft <= Lleft0 ||  Lright >= Lright0
            Lleft2 = L;
            L2 = L*1.5;
        end
        if Lright < Lright0 && Lleft > Lleft0
            L2 = 0.5*(Lright+ Lleft);
        end
        if Lleft > Lleft0 && L > Lleft
            Lleft2 = L;
        end
    end
    disp('Change Length Size');
else
    hasFind = 1;
end

end