function [Auxn,Auxs,Auyn,Auys,Axn,Axs,Ayn,Ays,Axyn,Axys]= CalcCoefficientMatrix_H(isTip,CurEle,Length,sinbeta,cosbeta,G,miu,x,y)
para = -1/4/pi/(1-miu);
% source element
a2 = Length(2)/2;

%Calculate angle of discontinuities
% relative locations
vecF = [cosbeta, sinbeta];%[CurEle(2,3)-CurEle(2,1),CurEle(2,4)-CurEle(2,2)]/(a2*2);
vecn = [-sinbeta,cosbeta];

% relative coordinates to the fracture
relR = [x-(CurEle(2,1)+CurEle(2,3))/2,y-(CurEle(2,2)+CurEle(2,4))/2];
relx = dot(relR,vecF);
rely = dot(relR,vecn);
if  abs(rely) < 1e-16
    rely = 0;
end
if  abs(relx) < 1e-16
    relx = 0;
end
%  Relative locations
theta1 = atan(rely/(relx-a2));
theta2 = atan(rely/(relx+a2));
r1 = sqrt((relx-a2)^2+rely^2);
r2 = sqrt((relx+a2)^2+rely^2);

% the angles are not right deal with
if relx == -a2 && rely > 0
    theta1 = theta1 + pi;
end
if relx == -a2 && rely < 0
    theta1 = theta1 - pi;
end
if abs(relx) < a2 && rely > 0
    theta1 = atan(rely/(relx-a2))+pi;
end
if abs(relx) < a2 && rely < 0
%     rely = -rely;
%     relx = -relx;
%     theta1 = atan(rely/(relx-a2))+pi;
%     theta2 = atan(rely/(relx+a2));
%     cosbeta = -cosbeta;
%     sinbeta = -sinbeta;
    theta1 = -pi + atan(rely/(relx-a2));
end
if relx < -a2 && rely > 0
    theta1 = pi + atan(rely/(relx-a2));
    theta2 = pi + atan(rely/(relx+a2));
end
if relx > a2 && rely < 0
    theta1 =  atan(rely/(relx-a2));
    theta2 =  atan(rely/(relx+a2));
end
if relx < -a2 && rely < 0
    theta1 = -pi + atan(rely/(relx-a2));
    theta2 = -pi + atan(rely/(relx+a2));
end

if isTip < 0.9
    Auxn0 = zeros(3,1);
    Auxs0 = Auxn0;
    Auyn0 = Auxn0;
    Auys0 = Auxn0;
    Axn0 = Auxn0;
    Axs0 = Auxn0;
    Ayn0 = Auxn0;
    Ays0 = Auxn0;
    Axyn0 = Auxn0;
    Axys0 = Auxn0;
    % This Element is not at tip of fractures
    a1 = Length(1)/2;
    a3 = Length(3)/2;
    
    % Derivatives
    I0x = CalcI(0,theta1,theta2,r1,r2,relx,rely,a2,2);
    I0y = CalcI(0,theta1,theta2,r1,r2,relx,rely,a2,3);
    I0xx = CalcI(0,theta1,theta2,r1,r2,relx,rely,a2,4);
    I0yy = CalcI(0,theta1,theta2,r1,r2,relx,rely,a2,5);
    I0xy = CalcI(0,theta1,theta2,r1,r2,relx,rely,a2,6);
    I0xyy = CalcI(0,theta1,theta2,r1,r2,relx,rely,a2,7);
    I0yyy = CalcI(0,theta1,theta2,r1,r2,relx,rely,a2,8);
    
    I1x = CalcI(1,theta1,theta2,r1,r2,relx,rely,a2,2);
    I1y = CalcI(1,theta1,theta2,r1,r2,relx,rely,a2,3);
    I1xx = CalcI(1,theta1,theta2,r1,r2,relx,rely,a2,4);
    I1yy = CalcI(1,theta1,theta2,r1,r2,relx,rely,a2,5);
    I1xy = CalcI(1,theta1,theta2,r1,r2,relx,rely,a2,6);
    I1xyy = CalcI(1,theta1,theta2,r1,r2,relx,rely,a2,7);
    I1yyy = CalcI(1,theta1,theta2,r1,r2,relx,rely,a2,8);
    
    I2x = CalcI(2,theta1,theta2,r1,r2,relx,rely,a2,2);
    I2y = CalcI(2,theta1,theta2,r1,r2,relx,rely,a2,3);
    I2xx = CalcI(2,theta1,theta2,r1,r2,relx,rely,a2,4);
    I2yy = CalcI(2,theta1,theta2,r1,r2,relx,rely,a2,5);
    I2xy = CalcI(2,theta1,theta2,r1,r2,relx,rely,a2,6);
    I2xyy = CalcI(2,theta1,theta2,r1,r2,relx,rely,a2,7);
    I2yyy = CalcI(2,theta1,theta2,r1,r2,relx,rely,a2,8);
    %DD high order element
    % F1 is the shape function of Left element
    F1y = 1/(a1+a2)/(a1+2*a2+a3)*(I2y - (a2+a3)*I1y);
    F2y = -1/(a1+a2)/(a2+a3)*(I2y + (a1-a3)*I1y) + I0y;
    F3y = 1/(a2+a3)/(a1+2*a2+a3)*(I2y + (a1+a2)*I1y);
    
    F1x = 1/(a1+a2)/(a1+2*a2+a3)*(I2x - (a2+a3)*I1x);
    F2x = -1/(a1+a2)/(a2+a3)*(I2x + (a1-a3)*I1x) + I0x;
    F3x = 1/(a2+a3)/(a1+2*a2+a3)*(I2x + (a1+a2)*I1x);
    
    F1xx = 1/(a1+a2)/(a1+2*a2+a3)*(I2xx - (a2+a3)*I1xx);
    F2xx = -1/(a1+a2)/(a2+a3)*(I2xx + (a1-a3)*I1xx) + I0xx;
    F3xx = 1/(a2+a3)/(a1+2*a2+a3)*(I2xx + (a1+a2)*I1xx);
    
    F1yy = 1/(a1+a2)/(a1+2*a2+a3)*(I2yy - (a2+a3)*I1yy);
    F2yy = -1/(a1+a2)/(a2+a3)*(I2yy + (a1-a3)*I1yy) + I0yy;
    F3yy = 1/(a2+a3)/(a1+2*a2+a3)*(I2yy + (a1+a2)*I1yy);
    
    F1xy = 1/(a1+a2)/(a1+2*a2+a3)*(I2xy - (a2+a3)*I1xy);
    F2xy = -1/(a1+a2)/(a2+a3)*(I2xy + (a1-a3)*I1xy) + I0xy;
    F3xy = 1/(a2+a3)/(a1+2*a2+a3)*(I2xy + (a1+a2)*I1xy);
    
    F1xyy = 1/(a1+a2)/(a1+2*a2+a3)*(I2xyy - (a2+a3)*I1xyy);
    F2xyy = -1/(a1+a2)/(a2+a3)*(I2xyy + (a1-a3)*I1xyy) + I0xyy;
    F3xyy = 1/(a2+a3)/(a1+2*a2+a3)*(I2xyy + (a1+a2)*I1xyy);
    
    F1yyy = 1/(a1+a2)/(a1+2*a2+a3)*(I2yyy - (a2+a3)*I1yyy);
    F2yyy = -1/(a1+a2)/(a2+a3)*(I2yyy + (a1-a3)*I1yyy) + I0yyy;
    F3yyy = 1/(a2+a3)/(a1+2*a2+a3)*(I2yyy + (a1+a2)*I1yyy);
    
    % stress and displacements in discontinuity coordinates
    Auxs0(1) = para*(2*(1-miu)*F1y - rely*F1xx);  %Correspond to DS(DX) Left elements
    Auxs0(2) = para*(2*(1-miu)*F2y - rely*F2xx);  %Correspond to DS(DX) source elements
    Auxs0(3) = para*(2*(1-miu)*F3y - rely*F3xx);  %Correspond to DS(DX) right elements
    
    Auxn0(1) = para*(-(1-2*miu)*F1x-rely*F1xy);   % Correspond to Dn(Dy)
    Auxn0(2) = para*(-(1-2*miu)*F2x-rely*F2xy);   % Correspond to Dn(Dy)
    Auxn0(3) = para*(-(1-2*miu)*F3x-rely*F3xy);   % Correspond to Dn(Dy)
    
    Auys0(1) = para*((1-2*miu)*F1x- rely*F1xy);
    Auys0(2) = para*((1-2*miu)*F2x- rely*F2xy);
    Auys0(3) = para*((1-2*miu)*F3x- rely*F3xy);
    
    Auyn0(1) = para*(2*(1-miu)*F1y-rely*F1yy);
    Auyn0(2) = para*(2*(1-miu)*F2y-rely*F2yy);
    Auyn0(3) = para*(2*(1-miu)*F3y-rely*F3yy);
    
    Axs0(1) = 2*G*para*(2*F1xy+rely*F1xyy);
    Axs0(2) = 2*G*para*(2*F2xy+rely*F2xyy);
    Axs0(3) = 2*G*para*(2*F3xy+rely*F3xyy);
    
    Axn0(1) = 2*G*para*(F1yy+rely*F1yyy);
    Axn0(2) = 2*G*para*(F2yy+rely*F2yyy);
    Axn0(3) = 2*G*para*(F3yy+rely*F3yyy);
    
    Ays0(1) = 2*G*para*(-rely*F1xyy) ;
    Ays0(2) = 2*G*para*(-rely*F2xyy) ;
    Ays0(3) = 2*G*para*(-rely*F3xyy) ;
    
    Ayn0(1) = 2*G*para*(F1yy-rely*F1yyy);
    Ayn0(2) = 2*G*para*(F2yy-rely*F2yyy);
    Ayn0(3) = 2*G*para*(F3yy-rely*F3yyy);
    
    Axys0(1) = 2*G*para*(F1yy+rely*F1yyy);
    Axys0(2) = 2*G*para*(F2yy+rely*F2yyy);
    Axys0(3) = 2*G*para*(F3yy+rely*F3yyy);
    
    Axyn0(1) = 2*G*para*(-rely*F1xyy);
    Axyn0(2) = 2*G*para*(-rely*F2xyy);
    Axyn0(3) = 2*G*para*(-rely*F3xyy);
    
    % $$$$$$ Coordinate transfer $$$$$
   % beta  = - beta;
    % Coordinate transformation
    % Displacement
    Auxs = Auxs0 * cosbeta - Auys0 * sinbeta;
    Auxn = Auxn0 * cosbeta - Auyn0 * sinbeta;
    Auys = Auxs0 * sinbeta + Auys0 * cosbeta;
    Auyn = Auxn0 * sinbeta + Auyn0 * cosbeta;
    % Stress
    %Crouch Page 24
    Axs = Axs0 * cosbeta* cosbeta - 2 * Axys0 * sinbeta * cosbeta + Ays0 * sinbeta * sinbeta;
    Axn = Axn0 * cosbeta* cosbeta - 2 * Axyn0 * sinbeta * cosbeta + Ayn0 * sinbeta * sinbeta;
    Axys = (Axs0 - Ays0) * sinbeta * cosbeta + Axys0 * (cosbeta* cosbeta-sinbeta * sinbeta);
    Axyn = (Axn0 - Ayn0) * sinbeta * cosbeta + Axyn0 * (cosbeta* cosbeta-sinbeta * sinbeta);
    Ays = Axs0 * sinbeta * sinbeta + 2 * Axys0 * sinbeta * cosbeta + Ays0 *  cosbeta* cosbeta;
    Ayn = Axn0 * sinbeta * sinbeta + 2 * Axyn0 * sinbeta * cosbeta + Ayn0 *  cosbeta* cosbeta;
end
end

