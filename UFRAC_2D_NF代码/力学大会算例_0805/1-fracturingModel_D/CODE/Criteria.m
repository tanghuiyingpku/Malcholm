function [initAngle,rc] = Criteria(angleHF,Sx0,Sy0,Sxy0,beta,KI,KII,To)
sinb = sind(angleHF);
cosb = cosd(angleHF);
Sy =  (Sx0*sinb^2-2*Sxy0*sinb*cosb + Sy0*cosb^2);
Sxy = (Sx0 - Sy0)*sinb*cosb + Sxy0*(cosb^2-sinb^2);
Sx = (Sx0*cosb^2+2*Sxy0*sinb*cosb + Sy0*sinb^2);

alpha = KII/KI;
%beta = 60;
a1 = (cosd(beta/2) - alpha*sind(beta/2))^2;
a2 = 2 *(Sx/2+Sy/2 - To)*(cosd(beta/2) - alpha*sind(beta/2));
a3 = (Sx/2+Sy/2-To)^2;
a4 = (Sx/2 - Sy/2)^2;
a5 = sind(beta/2)*cosd(beta/2)*sind(beta*3/2)+alpha*(sind(beta/2)+sind(beta/2)*cosd(beta/2)*cosd(beta*3/2));
a6 = sind(beta/2)*cosd(beta/2)*cosd(beta*3/2) + alpha*cosd(beta/2)*(1-sind(beta/2)*sind(beta*3/2));
A = a1 - a5^2 - a6^2;
B = a2 + (Sx- Sy)*a5 - 2*Sxy*a6;
C = a3 - a4 - Sxy^2;
k = (-B -sqrt(B^2 - 4*A*C))/2/A;
rc = (KI/k)^2/2/pi;
%
Tau = k*cosd(beta/2)*sind(beta/2)*cosd(beta*3/2) + k*alpha*cosd(beta/2)*(1-sind(beta/2)*sind(beta*3/2));
sigmaX = Sx + k*cosd(beta/2)*(1-sind(beta/2)*sind(beta*3/2)) + k*alpha*sind(beta/2)*(-2-cosd(beta/2)*cosd(beta*3/2));
sigmaY = Sy + k*cosd(beta/2)*(1+sind(beta/2)*sind(beta*3/2)) + k*alpha*sind(beta/2)*cosd(beta/2)*cosd(beta*3/2);
% Principle direction has two 
initAngle1 =0.5* atand(2*Tau/(sigmaX-sigmaY));
initAngle2 = (atand(2*Tau/(sigmaX-sigmaY))-180)/2;

sinb = sind(initAngle1);
cosb = cosd(initAngle1);
S1 =  (sigmaX*sinb^2-2*Tau*sinb*cosb + sigmaY*cosb^2);
sinb = sind(initAngle2);
cosb = cosd(initAngle2);
S2 =  (sigmaX*sinb^2-2*Tau*sinb*cosb + sigmaY*cosb^2);
if S1 > S2
    initAngle = initAngle1;
else
    initAngle = initAngle2;
end


end