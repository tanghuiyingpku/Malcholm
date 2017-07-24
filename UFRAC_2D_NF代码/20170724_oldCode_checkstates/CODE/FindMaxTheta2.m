function FindMaxTheta2(Kmax,Kmin)
Kic = 0;
a = -90;
b = 90;
coeffv = 0:2:5;
betav = 0:30:180;
%alpha KII/KI
angle = zeros(length(coeffv),length(betav));
for i = 1: length(coeffv)
    for j = 1:length(betav)
        coeff = coeffv(i);
        beta = betav(j);
        for theta0 = a : 0.1: b
            k1 = 0.5*cosd(theta0/2)*((1+cosd(theta0))-3*coeff*sind(theta0));
            k2 = 0.5*cosd(theta0/2)*(sind(theta0)+coeff*(3*cosd(theta0)-1));
            % Angle between Kmax and K
            alpha = beta - theta0 ;
            if abs(alpha) < 1e-8
                Kalpha = Kmax;
            else
                kk = cotd(alpha);
                Kalpha = sqrt(1+kk^2)*sqrt(1/(kk^2/Kmax^2+1/Kmin^2));
            end
            Ki = (k1^2+k2^2)/Kalpha^2;
            if Ki > Kic
                angle(i,j) = theta0;
            end
        end
    end
end
save ANGLE.txt -ascii angle;
end