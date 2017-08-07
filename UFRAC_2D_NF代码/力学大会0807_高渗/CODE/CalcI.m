function II = CalcI(Type,theta1,theta2,r1,r2,relx,rely,a,i)
%I0
if abs(Type)<1e-6
    II = CalcI0(r1,r2,theta1,theta2,relx,rely,a,i);
    return;
end
if abs(Type-1)<1e-6
    II = CalcI1(r1,r2,theta1,theta2,relx,rely,a,i);
    return;
end
if abs(Type-2)<1e-6
    II = CalcI2(r1,r2,theta1,theta2,relx,rely,a,i);
    return;
end
if Type-3 > -1e-6  % two possibilities : left tip and right tip
    II = CalcIC(relx,rely,a,i,Type);
    return;

end
end
function II = CalcI0(r1,r2,theta1,theta2,relx,rely,a,i)
if i == 1
    %Io
    II = rely*(theta1-theta2)-(relx-a)*log(r1)+(relx+a)*log(r2)-2*a;
    return;
end
if i == 2
    %Iox
    II = -log(r1/r2);
    return;
end
if i == 3
    %Ioy
    II = theta1 - theta2;
    return;
end
if i == 4
    %Ioxx
    II = -((relx-a)/r1^2-(relx+a)/r2^2);
    return;
end
if i == 5
    %Ioyy
    II = ((relx-a)/r1^2-(relx+a)/r2^2);
    return;
end
if i == 6
    %Ioxy
    II = -(rely/r1^2-rely/r2^2);
    return;
end
if i ==7
    %Ioxyy
    II = -((relx-a)^2-rely^2)/r1^4+((relx+a)^2-rely^2)/r2^4;
    return;
end
if i == 8
    %Ioyyy
    II = -(2*(relx-a)*rely)/r1^4+(2*(relx+a)*rely)/r2^4;
    return;
end
end
function II = CalcI1(r1,r2,theta1,theta2,relx,rely,a,i)
if i == 1
    %I1
    II = relx*rely*(theta1-theta2) + 0.5 * ( rely^2 - relx^2 + a^2)*(log(r1)-log(r2))-relx*a;
    return;
end
if i == 2
    %I1x
    II = rely*(theta1 - theta2) - relx*log(r1/r2) - 2*a;
    return;
end
if i == 3
    %I1y
    II = relx*(theta1 - theta2) + rely*log(r1/r2);
    return;
end
if i == 4
    %I1xx
    II = -(log(r1/r2) + a*(relx-a)/r1^2 + a*(relx+a)/r2^2);
    return;
end
if i == 5
    %I1yy
    II = (log(r1/r2) + a*(relx-a)/r1^2 + a*(relx+a)/r2^2);
    return;
end
if i == 6
    %I1xy
    %  II=0.3478548*f1xy(0.8611363,relx,rely,a) + 0.3478548*f1xy(-0.8611363,relx,rely,a) + 0.6521452*f1xy(0.3399810,relx,rely,a)+0.6521452*f1xy(-0.3399810,relx,rely,a);
    % II = 0.4679139*f1xy(0.2386192,relx,rely,a)+
    % 0.4679139*f1xy(-0.2386192,relx,rely,a)+ 0.3607616*f1xy(0.6612094,relx,rely,a)+0.3607616*f1xy(-0.6612094,relx,rely,a)+ 0.1713245*f1xy(0.9324695,relx,rely,a)+ 0.1713245*f1xy(-0.9324695,relx,rely,a);
    II = (theta1 - theta2) -2*a*rely*(a^2+relx^2+rely^2)/((a^2-relx^2)^2+2*(a^2+relx^2)*rely^2+rely^4);
    return;
end
if i ==7
    %I1xyy
    II = relx / r1^2 - relx/r2^2 - 2*a*(relx - a)^2/r1^4 - 2*a*(relx + a)^2/r2^4;
    return;
end
if i == 8
    %I1yyy
    II = rely / r1^2 - rely/r2^2 - 2*a*(relx - a)*rely/r1^4 - 2*a*(relx + a)*rely/r2^4;
    return;
end
end

function II = CalcI2(r1,r2,theta1,theta2,x,y,a,i)
if i == 1
    %I2
    II = (y/3)*(3*x^2 - y^2)*(theta1-theta2) + 1/3*(3*x*y^2 - x^3 + a^3)*log(r1) - 1/3*(3*x*y^2 - x^3 -a^3)*log(r2) - 2/3*a*(x^2 - y^2 + a^2/3);
    return;
end
if i == 2
    %I2x
    II = 2*x*y*(theta1 - theta2) +(y^2 - x^2)*log(r1/r2) - 2*a*x;
    return;
end
if i == 3
    %I2y
    II = (x^2 - y^2)*(theta1 - theta2) + 2*x*y*log(r1/r2) + 2*a*y;
    return;
end
if i == 4
    %I2xx
    II = 2*y*(theta1 - theta2) - 2*x*log(r1/r2) - a^2*(x-a)/r1^2 + a^2*(x+a)/r2^2 - 4*a;
    return;
end
if i == 5
    %I2yy
    II = -(2*y*(theta1 - theta2) - 2*x*log(r1/r2) - a^2*(x-a)/r1^2 + a^2*(x+a)/r2^2 - 4*a);
    return;
end
if i == 6
    %I2xy
    II = 2*x*(theta1 -theta2) + 2*y*log(r1/r2) - a^2*y/r1^2 + a^2*y/r2^2;
    return;
end
if i ==7
    %I2xyy
    II = 2*log(r1/r2) + a*(2*x - a)/r1^2 + a*(2*x+a)/r2^2 -2*a^2*(x-a)^2/r1^4 + 2*a^2*(x+a)^2/r2^4;
    return;
end
if i == 8
    %I2yyy
    II = -2*(theta1 - theta2) + 2*a*y/r1^2 + 2*a*y/r2^2 -2*a^2*(x-a)*y/r1^4 + 2*a^2*(x+a)*y/r2^4;
    return;
end
end
function II = CalcIC(x,y,a,i,Type)
% In this coordinate system start from Tip top

if i == 2
    %Icx
    II = CalcB(2,x,y,a);
    if Type > 3.1
        II = II * (-1);
    end
    %II = x*(A1(1) - A1(2)) - (A2(1) - A2(2));
    %     II1 = quadgk(@(e)sqrt(e).*(x-e)./((x-e).^2+y^2),0,a0-1e-8);
    %     II2 = quadgk(@(e)sqrt(e).*(x-e)./((x-e).^2+y^2),a0+1e-8,2*a0);
    %     II = II1+II2;
    return;
end
if i == 3
    %Icy
    II = CalcB(3,x,y,a);
    %     II1 = quadgk(@(e)sqrt(e)*y./((x-e).^2+y^2),0,a0-1e-8);
    %     II2 = quadgk(@(e)sqrt(e)*y./((x-e).^2+y^2),a0+1e-8,2*a0);
    %     II = II1+II2;
    return;
    
end
if i == 4
    %Icxx
    II = CalcB(5,x,y,a);
    %     II1 = quadgk(@(e)(-sqrt(e)*y./((x-e).^2+y^2)+2*sqrt(e)*y^2./((x-e).^2+y^2).^2),0,a0-1e-8);
    %     II2 = quadgk(@(e)(-sqrt(e)*y./((x-e).^2+y^2)+2*sqrt(e)*y^2./((x-e).^2+y^2).^2),a0+1e-8,2*a0);
    %     II = II1+II2;
    return;
    
end
if i == 5
    %Icyy
    II = -CalcB(5,x,y,a);
    %     II1 = quadgk(@(e)(sqrt(e)./((x-e).^2+y^2)-2*sqrt(e)*y^2./((x-e).^2+y^2).^2),0,a0-1e-8);
    %     II2 = quadgk(@(e)(sqrt(e)./((x-e).^2+y^2)-2*sqrt(e)*y^2./((x-e).^2+y^2).^2),a0+1e-8,2*a0);
    %     II = II1+II2;
    return;
end
if i == 6
    %Icxy
    II = CalcB(4,x,y,a);
    if Type > 3.1
        II = II * (-1);
    end
    %     II1 = quadgk(@(e)(-2*sqrt(e)*y.*(x-e)./((x-e).^2+y^2).^2),0,a0-1e-8);
    %     II2 = quadgk(@(e)(-2*sqrt(e)*y.*(x-e)./((x-e).^2+y^2).^2),a0+1e-8,2*a0);
    %     II = II1+II2;
    return;
    
end
if i ==7
    %Icxyy
    II = CalcB(6,x,y,a);
    if Type > 3.1
        II = II * (-1);
    end
    %     II1 = quadgk(@(e)(-2*sqrt(e).*(x-e)./((x-e).^2+y^2).^2+8*sqrt(e)*y^2.*(x-e)./((x-e).^2+y^2).^3),0,a0-1e-8);
    %     II2 = quadgk(@(e)(-2*sqrt(e).*(x-e)./((x-e).^2+y^2).^2+8*sqrt(e)*y^2.*(x-e)./((x-e).^2+y^2).^3),a0+1e-8,2*a0);
    %     II = II1+II2;
    return;
    
end
if i == 8
    %Icyyy
    II = CalcB(7,x,y,a);
    %     II1 = quadgk(@(e)(-6*sqrt(e)*y./((x-e).^2+y^2).^2+8*sqrt(e)*y^3./((x-e).^2+y^2).^3),0,a0-1e-8);
    %     II2 = quadgk(@(e)(-6*sqrt(e)*y./((x-e).^2+y^2).^2+8*sqrt(e)*y^3./((x-e).^2+y^2).^3),a0+1e-8,2*a0);
    %     II = II1+II2;
    return;
    
end
end
function B = CalcB(i,x,y,a)
eps = 1e-8;
if abs(y) < eps
    if abs(i-4) < eps
        B = 0;
        return;
    end
    if abs(i-7) < eps
        B = 0;
        return;
    end
    if x > -a
        if abs(i-2) < eps
            B = -2*sqrt(2)+sqrt((x+a)/a)*log(abs((sqrt(x+a)+sqrt(2*a))/(sqrt(x+a)-sqrt(2*a))));
            return;
        end
        if abs(i-5) < eps
          %THY  B = (-1)*(1/2/sqrt(a*(a+x))*log(abs((sqrt(x+a)+sqrt(2*a))/(sqrt(x+a)-sqrt(2*a))))-sqrt(2)/(x-a));
          B = (1)*(1/2/sqrt(a*(a+x))*log(abs((sqrt(x+a)+sqrt(2*a))/(sqrt(x+a)-sqrt(2*a))))-sqrt(2)/(x-a));
          return;
        end
        if abs(i-6) < eps
            %THY B = (1)*(sqrt(2)/2/(x^2-a^2)+1/4/sqrt(a)/(x+a)^1.5*log(abs((sqrt(x+a)+sqrt(2*a))/(sqrt(x+a)-sqrt(2*a))))-sqrt(2)/(x-a)^2);
            B = (1)*(sqrt(2)/2/(x^2-a^2)+1/4/sqrt(a)/(x+a)^1.5*log(abs((sqrt(x+a)+sqrt(2*a))/(sqrt(x+a)-sqrt(2*a))))-sqrt(2)/(x-a)^2);
            return;
        end
    else
        r = abs(x)-a;
        if abs(i-2) < eps
            B = -2*sqrt(2)+2*sqrt(r/a)*atan(sqrt(2*a/r));
            return;
        end
        if abs(i-5) < eps
            B = sqrt(2)/(r+2*a)-1/sqrt(a*r)*atan(2*a/r);
            return;
        end
        if abs(i-6) < eps
            B = sqrt(2)/2/r/(r+2*a) + 1/2/sqrt(a)/r^1.5*atan(2*a/r)-sqrt(2)/(r+2*a)^2;
            return;
        end
    end
    if abs(i-3) < eps
        if abs(x) > a
            B=0;
        else
            %THY 2014/10/13
            if y > 0
                B = pi;
            else
                B = -pi;
            end
%             if y > 0
%                 B = -pi;
%             else
%                 B = pi;
%             end
        end
        return;
    end
else
    % B = 0.3478548*f(i,0.8611363,x,y,a) + 0.3478548*f(i,-0.8611363,x,y,a) + 0.6521452*f(i,0.3399810,x,y,a)+0.6521452*f(i,-0.3399810,x,y,a);
    B = 0.4679139*f(i,0.2386192,x,y,a)+ 0.4679139*f(i,-0.2386192,x,y,a)+ 0.3607616*f(i,0.6612094,x,y,a)+...
        0.3607616*f(i,-0.6612094,x,y,a)+ 0.1713245*f(i,0.9324695,x,y,a)+ 0.1713245*f(i,-0.9324695,x,y,a);
end
end
function y=f(i,t,x,y,a)
if abs(i-2) < 1e-6
    y = a*(x-a*t)/((x-a*t)^2+y^2)*(1+t)^0.5;
    return;
end
if abs(i-3) < 1e-6
    y = a*y/((x-a*t)^2+y^2)*(1+t)^0.5;
    return;
end
if abs(i-4) < 1e-6
    y = -a*2*y*(x-a*t)/((x-a*t)^2+y^2)^2*(1+t)^0.5;
    return;
end
if abs(i-5) < 1e-6
    y = -a*((x-a*t)^2-y^2)/((x-a*t)^2+y^2)^2*(1+t)^0.5;
    return;
end
if abs(i-6) < 1e-6
    y = -2*a*((x-a*t)^3/((x-a*t)^2+y^2)^3-3*(x-a*t)*y^2/((x-a*t)^2+y^2)^3)*(1+t)^0.5;
    return;
end
if abs(i-7) < 1e-6
    y = -2*a*y*(3*(x-a*t)^2/((x-a*t)^2+y^2)^3-y^2/((x-a*t)^2+y^2)^3)*(1+t)^0.5;
    return;
end
end
