function tempCode
r = stepL * 0.01;
        sinb = AllElenew(32,5);
        cosb = AllElenew(32,6);
        xy0 = [AllEle(32,1),AllEle(32,2)];
        angle = 150:1:210;
        Sxi = zeros(length(angle),1);
        Syi = zeros(length(angle),1);
        Sxyi = zeros(length(angle),1);
        Sn = zeros(length(angle),1);
        Tt = zeros(length(angle),1);
        shh = Sn;
        sHH = Sn;
        for ii = 1 : length(angle)
            ca = angle(ii);
            cosbi = cosb*cosd(ca) - sinb*sind(ca);
            sinbi = sinb*cosd(ca)+cosb*sind(ca);
            xyi = xy0 + r*[cosbi,sinbi];
            [Sxi(ii),Syi(ii),Sxyi(ii),~,~] = CalcPointStress_C(Mat,DD_new,xyi(1),xyi(2),AllElenew,ConnList);
            Sxi(ii) = Sxi(ii) + Sxx;
            Syi(ii) = Syi(ii) + Syy;
            Sn(ii) = Sxi(ii)*sinbi^2-2*Sxyi(ii)*sinbi*cosbi + Syi(ii)*cosbi^2;
            Tt(ii)= -(Sxi(ii) - Syi(ii))*sinbi*cosbi+ Sxyi(ii)*(cosbi^2 - sinbi^2);
            [shh(ii),sHH(ii)] =  PrinStress2(Sxi(ii),Syi(ii),Sxyi(ii));
        end
        
    [a,b] = max(Sn);
    angle(b)
        
        
        %%
         r = stepL * 0.05;
    sinb = 1;
    cosb = 0;
    xy0 = [-2,1];
    angle = 6:0.05:9;
    Sxi = zeros(length(angle),1);
    Syi = zeros(length(angle),1);
    Sxyi = zeros(length(angle),1);
    Sn = zeros(length(angle),1);
    Tt = zeros(length(angle),1);
    shh = Sn;
    sHH = Sn;
    figure;
    plot(xy0(1),xy0(2),'r*');
    for ii = 1 : length(angle)
        ca = angle(ii);
        cosbi = cosb*cosd(ca) - sinb*sind(ca);
        sinbi = sinb*cosd(ca)+cosb*sind(ca);
        xyi = xy0 + r*[cosbi,sinbi];
        plot(xyi(1),xyi(2),'k*');
        hold on;
        [Sxi(ii),Syi(ii),Sxyi(ii),~,~] = CalcPointStress_C(Mat,DD_new,xyi(1),xyi(2),AllElenew,ConnList);
        Sxi(ii) = Sxi(ii) + Sxx;
        Syi(ii) = Syi(ii) + Syy;
        Sn(ii) = Sxi(ii)*sinbi^2 - 2*Sxyi(ii)*sinbi*cosbi + Syi(ii)*cosbi^2;
%         Tt(ii)= (Sxi(ii) - Syi(ii))*sinbi*cosbi+ Sxyi(ii)*(cosbi^2 - sinbi^2);
%         shh(ii) = Sxi(ii)*cosbi^2-2*Sxyi(ii)*sinbi*cosbi + Syi(ii)*sinbi^2;
    end
    % Stress infront of tip
    [a,b] = max(Sn);
    angle(b)
end