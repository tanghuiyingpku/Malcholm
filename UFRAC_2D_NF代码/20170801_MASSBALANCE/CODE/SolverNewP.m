function newP = SolverNewP(nFracture,Fractures,nAllEle,AllEle,ConnList,Mat,Pres,DD)
global GrowthN
global stepL;
newP = zeros(GrowthN,1);
Pres = Pres*1e6;
BC   = zeros(2*nAllEle,1);
InsituS = zeros(nAllEle*2,1);
% Project the in-situ stresses to fractures;
Sxx = Mat.Sxx*1e6;
Syy = Mat.Syy*1e6;
Sxy = Mat.Sxy*1e6;
for i = 1 : nAllEle
    %alpha = theta + 90;
    sinalp = AllEle(i,6);  %cosbet(i);
    cosalp = -AllEle(i,5); %-sinbet(i);
    InsituS(i) = (Sxx - Syy) * cosalp * sinalp - Sxy*(cosalp^2 - sinalp^2);
    InsituS(i + nAllEle) = Sxx * cosalp^2 + Syy * sinalp^2+2*Sxy*sinalp*cosalp;
end
BC_insitu = InsituS;
for ii = 1 : nAllEle
    BC(ii) = BC(ii)-BC_insitu(ii);
    BC(nAllEle+ii) =BC(nAllEle+ii) -BC_insitu(nAllEle+ii);
    % Derivative of Pressure
end
BC(nAllEle+1:nAllEle*2-GrowthN) = BC(nAllEle+1:nAllEle*2-GrowthN) - Pres;
[~,CM]  = BuildCoefMatix_H2(Mat,nFracture,Fractures,nAllEle,AllEle,ConnList,1); 
save CM0.txt -ascii CM;
% The unknowns change to be Pressure
A = inv(CM);
B = zeros(nAllEle*2,1);
B(1:nAllEle) = [DD(1:nAllEle-GrowthN);zeros(GrowthN,1)];
%Square Decrease
openTip = DD((nAllEle-GrowthN)*2);
Ltip = AllEle(nAllEle-GrowthN,7);
x = 0:GrowthN-1;
xx = [x'*stepL;stepL*(GrowthN-0.5)+Ltip/2];
alpha = openTip/sqrt(xx(GrowthN+1));
openD = alpha * sqrt(abs(xx(GrowthN:-1:1)));


B(nAllEle+1:2*nAllEle) =[DD((nAllEle-GrowthN)+1:(nAllEle-GrowthN)*2);openD];% [DD((nAllEle-GrowthN)+1:(nAllEle-GrowthN)*2);zeros(GrowthN,1)];
%New Pressure Distribution
XX = A\B;
for i = 1 : GrowthN
    newP(i) = -XX(nAllEle*2-GrowthN+i) - InsituS(nAllEle*2-GrowthN+i);
    newP(i) = newP(i)/1e6;
end
end
