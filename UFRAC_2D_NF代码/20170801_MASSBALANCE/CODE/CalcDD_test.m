function DD = CalcDD_test(Mat,nEle,Pres,AllEle)
Pres = Pres*1e6;
% Project the in-situ stresses to fractures;
Sxx = Mat.Sxx*1e6;
Syy = Mat.Syy*1e6;
Sxy = Mat.Sxy*1e6;

InsituS = zeros(nEle*2,1);
BC = zeros(nEle*2,1);
for i = 1 : nEle
    %alpha = theta + 90;
    sinalp = AllEle(i,6);  %cosbet(i);
    cosalp = -AllEle(i,5); %-sinbet(i);
    InsituS(i) = (Sxx - Syy) * cosalp * sinalp - Sxy*(cosalp^2 - sinalp^2);
    InsituS(i + nEle) = Sxx * cosalp^2 + Syy * sinalp^2+2*Sxy*sinalp*cosalp;
end
BC_insitu = InsituS;

for ii = 1 : nEle
    BC(ii) = BC(ii)-BC_insitu(ii);
    BC(nEle+ii) =BC(nEle+ii) -BC_insitu(nEle+ii) - Pres(ii);
end
disp('           ');
[extraBC,CM]  = BuildCoefMatix_Constant_test(nEle,Mat,AllEle);
%Unkown Displacements
%Solve Linear System for Displacements
disp('           ');
disp('Matrix Solving...');
disp('           ');
BC = BC + extraBC;
DD = CM\BC;

end