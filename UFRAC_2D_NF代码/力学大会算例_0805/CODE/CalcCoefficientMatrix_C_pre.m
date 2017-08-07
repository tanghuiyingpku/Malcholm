function [Auxn,Auxs,Auyn,Auys,Axn,Axs,Ayn,Ays,Axyn,Axys]= CalcCoefficientMatrix_C_pre(isTip,G,miu,FracLoc,x,y,sinbeta,cosbeta,length)
% relative locations
% vecF = [FracLoc(3)-FracLoc(1),FracLoc(4)-FracLoc(2)]/length;
% vecn = [-vecF(2),vecF(1)];
vecF = [cosbeta,sinbeta];%[CurEle(2,3)-CurEle(2,1),CurEle(2,4)-CurEle(2,2)]/(a2*2);
vecn = [-sinbeta,cosbeta];
%Calculate angle of discontinuities

% relative coordinates to the fracture
relR = [x-(FracLoc(1)+FracLoc(3))/2,y-(FracLoc(2)+FracLoc(4))/2];
relx = dot(relR,vecF);
rely = dot(relR,vecn);
a = length/2;

%  Relative locations
theta1 = atan(rely/(relx-a));
theta2 = atan(rely/(relx+a));
r1 = sqrt((relx-a)^2+rely^2);
r2 = sqrt((relx+a)^2+rely^2);
% 
if abs(relx) < a && rely > 0
    theta1 = atan(rely/(relx-a))+pi;
    theta2 = atan(rely/(relx+a));
end
if abs(relx) < a && rely < 0
    theta1 = pi + atan(rely/(relx-a));
    theta2 = 2*pi+atan(rely/(relx+a));
end
if relx >= a && rely > 0
    theta1 = atan(rely/(relx-a));
    theta2 = atan(rely/(relx+a));
end
if relx <= -a && rely > 0
    theta1 = pi + atan(rely/(relx-a));
    theta2 = pi + atan(rely/(relx+a));
end
if relx >= a && rely < 0
    theta1 =  atan(rely/(relx-a));
    theta2 =  atan(rely/(relx+a));
end
if relx <= -a && rely < 0
    theta1 = -pi + atan(rely/(relx-a));
    theta2 = -pi + atan(rely/(relx+a));
end

if isTip < 0.1
    %CDD element
    Ix = CalcI(0,theta1,theta2,r1,r2,relx,rely,a,2);
    Iy = CalcI(0,theta1,theta2,r1,r2,relx,rely,a,3);
    Ixx = CalcI(0,theta1,theta2,r1,r2,relx,rely,a,4);
    Iyy = CalcI(0,theta1,theta2,r1,r2,relx,rely,a,5);
    Ixy = CalcI(0,theta1,theta2,r1,r2,relx,rely,a,6);
    Ixyy = CalcI(0,theta1,theta2,r1,r2,relx,rely,a,7);
    Iyyy = CalcI(0,theta1,theta2,r1,r2,relx,rely,a,8);
else
    % Do not rotate any more but change the sign
    if isTip > 1.1
        relx = - relx;
    end
    % specil Tip element square root type
    Ix = CalcI(isTip+2,theta1,theta2,r1,r2,relx,rely,a,2);
    Iy = CalcI(isTip+2,theta1,theta2,r1,r2,relx,rely,a,3);
    Ixx = CalcI(isTip+2,theta1,theta2,r1,r2,relx,rely,a,4);
    Iyy = CalcI(isTip+2,theta1,theta2,r1,r2,relx,rely,a,5);
    Ixy = CalcI(isTip+2,theta1,theta2,r1,r2,relx,rely,a,6);
    Ixyy = CalcI(isTip+2,theta1,theta2,r1,r2,relx,rely,a,7);
    Iyyy = CalcI(isTip+2,theta1,theta2,r1,r2,relx,rely,a,8);
end
para = -1/4/pi/(1-miu);

% stress and displacements in discontinuity coordinates
Auxs0 = para*(2*(1-miu)*Iy - rely*Ixx);  %Correspond to Ds(Dx)
Auxn0 = para*(-(1-2*miu)*Ix-rely*Ixy);   % Correspond to Dn(Dy)
Auys0 = para*((1-2*miu)*Ix- rely*Ixy);
Auyn0 = para*(2*(1-miu)*Iy-rely*Iyy);
Axs0 = 2*G*para*(2*Ixy+rely*Ixyy);
Axn0 = 2*G*para*(Iyy+rely*Iyyy);
Ays0 = 2*G*para*(-rely*Ixyy) ;
Ayn0 = 2*G*para*(Iyy-rely*Iyyy);
Axys0 = 2*G*para*(Iyy+rely*Iyyy);
Axyn0 = 2*G*para*(-rely*Ixyy);

% the element is rotated so the coordinate should also rotate
% if isTip > 1.1
%     sinbeta = -sinbeta;
%     cosbeta = -cosbeta;
% end
% Coordinate transformation
% Displacement
%$$$$$$$$$$$ check  &$$$$$$$$$$$$$$
%beta = -beta;
%$$$$$$$$$$ anti-clockwise circle?????

Auxs = Auxs0 * cosbeta - Auys0 * sinbeta;
Auxn = Auxn0 * cosbeta - Auyn0 * sinbeta;
Auys = Auxs0 * sinbeta + Auys0 * cosbeta;
Auyn = Auxn0 * sinbeta + Auyn0 * cosbeta;
% Stress
Axs = Axs0 * cosbeta* cosbeta - 2 * Axys0 * sinbeta * cosbeta + Ays0 * sinbeta * sinbeta;
Axn = Axn0 * cosbeta* cosbeta - 2 * Axyn0 * sinbeta * cosbeta + Ayn0 * sinbeta * sinbeta;
Axys = (Axs0 - Ays0) * sinbeta * cosbeta + Axys0 * (cosbeta* cosbeta-sinbeta * sinbeta);
Axyn = (Axn0 - Ayn0) * sinbeta * cosbeta + Axyn0 * (cosbeta* cosbeta-sinbeta * sinbeta);
Ays = Axs0 * sinbeta * sinbeta + 2 * Axys0 * sinbeta * cosbeta + Ays0 *  cosbeta* cosbeta;
Ayn = Axn0 * sinbeta * sinbeta + 2 * Axyn0 * sinbeta * cosbeta + Ayn0 *  cosbeta* cosbeta;
end

