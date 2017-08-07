function ds = M_C_stability(angleHF,Sx0,Sy0,Sxy0,beta,KI,KII,rc,Fang,So)
P = 0;
sinb = sind(angleHF);
cosb = cosd(angleHF);
Sy =  (Sx0*sinb^2-2*Sxy0*sinb*cosb + Sy0*cosb^2);
Sxy = (Sx0 - Sy0)*sinb*cosb + Sxy0*(cosb^2-sinb^2);
Sx = (Sx0*cosb^2+2*Sxy0*sinb*cosb + Sy0*sinb^2);

k1 = KI/sqrt(rc*2*pi);
k2 = KII/sqrt(rc*2*pi);
Tau = Sxy + k1*cosd(beta/2)*sind(beta/2)*cosd(beta*3/2) + k2*cosd(beta/2)*(1-sind(beta/2)*sind(beta*3/2));
sigmaX = Sx + k1*cosd(beta/2)*(1-sind(beta/2)*sind(beta*3/2)) + k2*sind(beta/2)*(-2-cosd(beta/2)*cosd(beta*3/2));
sigmaY = Sy + k1*cosd(beta/2)*(1+sind(beta/2)*sind(beta*3/2)) +k2*sind(beta/2)*cosd(beta/2)*cosd(beta*3/2);
Sn = sigmaX*sind(beta)^2 - 2 * Tau * sind(beta)*cosd(beta) + sigmaY*cosd(beta)^2;
Sn = Sn + P;
t = (sigmaY - sigmaX)*sind(beta)*cosd(beta) + Tau *((cosd(beta))^2 - (sind(beta))^2);
ds = - abs(t) - tand(Fang)*Sn + So;
end