function StressLine()
close all
ratio = 0:0.04:1;
Crit = zeros(1,length(ratio));
Syh = -50;
Fangv = 0.05:0.002:1.0;
iscross = zeros(length(Fangv),1);
So = 0;
L = 0.5;
Sxyh = 0;
KI = 1;
KII = 0.0;
d = 0.0 ;
alphav =[30,45,60,75,90];
color = {'C','M','B','R','K'};
li = L/2;
x0 = d;
y0 = 0;
N = 50;
temp = zeros(1,length(ratio));
rt = zeros(length(alphav),1);
for ia = 1 : length(alphav)
    alpha = alphav(ia);
    for is =  1: length(ratio)
        fprintf('Stress Contrast = %f\n',10.^ratio(is));
        for ic = 1 : length(Fangv)
            Sxh = 10^ratio(is)*Syh;
            Fang = atand(Fangv(ic));
            xx= linspace(-li*cosd(alpha),li*cosd(alpha),N)+x0;
            yy = linspace(-li*sind(alpha),li*sind(alpha),N);
            rv = zeros(N,1);
            sinb  = sind(alpha);
            cosb = cosd(alpha);
            sigmaN = zeros(N,1);
            Sx = zeros(N,1);
            Sy = Sx;
            Sxy = Sx;
            tang = zeros(N,1);
            tau = zeros(N,1);
            Sh = zeros(N,1);
            Shmax = zeros(N,1);
            for i = 1 : N
                xi = xx(i);
                yi = yy(i);
                vec0 = [1 0];
                veci = [xi+d,yi];
                cosdi = dot(vec0,veci)/sqrt(veci(1)^2+veci(2)^2);
                theta = acosd(cosdi);
                r = sqrt(xi^2+yi^2);
                if yi < -1e-16
                    theta = -theta;
                    rv(i) = -r;
                    
                else
                    theta = theta;
                    rv(i) = r;
                end
                
                Sxxi = 1/sqrt(2*pi*r)*(KI*cosd(theta/2)*(1-sind(theta/2)*sind(3*theta/2))-KII*sind(theta/2)*(2+cosd(theta/2)*cosd(theta*3/2)));
                Syyi = 1/sqrt(2*pi*r)*(KI*cosd(theta/2)*(1+sind(theta/2)*sind(3*theta/2))+KII*sind(theta/2)*cosd(theta/2)*cosd(theta*3/2));
                Sxyi =1/sqrt(2*pi*r)*(KI*sind(theta/2)*cosd(theta/2)*cosd(3*theta/2)+KII*cosd(theta/2)*(1-sind(theta/2)*sind(theta*3/2)));
                Sxxi = Sxxi + Sxh;
                Syyi = Syyi + Syh;
                Sxyi = Sxyi + Sxyh;
                Sy(i) = Syyi;
                Sx(i) = Sxxi;
                Sxy(i) = Sxyi;
                Shmax(i) = (Sxxi+Syyi)/2+sqrt((Sxxi-Syyi)^2/4+Sxyi^2);
                Sh(i) = (Sxxi+Syyi)/2-sqrt((Sxxi-Syyi)^2/4+Sxyi^2);
                sigmaN(i)=  (Sxxi*sinb^2-2*Sxyi*sinb*cosb + Syyi*cosb^2);
                tau(i)= (Sxxi - Syyi)*sinb*cosb + Sxyi*(cosb^2-sinb^2);
                tang(i)= (Sxxi*cosb^2+2*Sxyi*sinb*cosb + Syyi*sinb^2);
            end
            
            Tn = 0;
            [anlger,rcr]= Criteria(0,Sxh,Syh,Sxyh,alpha,KI,KII,Tn);
            %         plot(rcr,0,'k.','Markersize',15);
            [anglel,rcl]= Criteria(0,Sxh,Syh,Sxyh,alpha - 180,KI,KII,Tn);
            %         plot(-rcl,0,'g.','Markersize',15);
            rt(i) = max(rcr,rcl);
            dSr = M_C_stability(0,Sxh,Syh,Sxyh,alpha,KI,KII,rcr,Fang,So);
            %dSr > 0 stable
            dSl = M_C_stability(0,Sxh,Syh,Sxyh,alpha-180,KI,KII,rcl,Fang,So);
            hasfind = 0;
            Crit(is) = nan;
            if min(dSr,dSl) < 1e-16
                iscross(ic) = 0;
                % legend('NO CROSSING');
            else
                %Cross
                hasfind = 1;
                Crit(is) = tand(Fang);
                fprintf('Critical Friction  = %f\n',Crit(is));
                iscross(ic) = 1;
                % legend('Cossing');
            end
            if hasfind > 0.1
                break;
            end
        end
    end
    temp = [temp;Crit];
    if ia == 1
        semilogy(Crit,10.^ratio,'color',[0.0,0.8,0.2],'Linewidth',2.5);
    else
        semilogy(Crit,10.^ratio,[color{ia}],'Linewidth',2.5);
    end
    hold on;
end
rt
legend('15','30','45','60','75','90');
save Criteria_P0.txt -ascii temp;
% figure
% plot(10.^ratio,Crit,'m*');
% plot(rv,sigmaN,'b','Linewidth',2.5);
% hold on
% plot(rv,tau,'r.');
% plot(rv,tang,'c--','Linewidth',2.5);
% legend('Normal','Shear','Tangential');
% plot([250 250],[-30 10],'m--');
end
