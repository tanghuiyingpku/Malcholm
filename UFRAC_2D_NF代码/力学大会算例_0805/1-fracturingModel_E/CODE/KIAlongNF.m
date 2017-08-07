function  [K11,theta11,K22,theta22,K33,theta33] = KIAlongNF()
global  Index IndexInv nAct DD_global TipStates;
global  Mat  AllEle_global ConnList_global MaxEle  ;
K11=0;
K22=0;
K33=0;
theta11 =0;
theta22 = 0;
theta33 = 0;
nEle = nAct;
DD = zeros(nEle*2,1);
AllEle = zeros(nEle,11);
ConnList = zeros(nEle,8);
num = zeros(3,1);
for i = 1 : nAct
    DD(i) = DD_global(IndexInv(i));
    DD(i+nEle) = DD_global(MaxEle+IndexInv(i));
    AllEle(i,:)  = AllEle_global(IndexInv(i),:);
    ConnList(i,:) = ConnList_global(IndexInv(i),:);
    if TipStates(IndexInv(i)) == -11
        for kk = 1 : ConnList_global(IndexInv(i),2)
            nb = ConnList_global(IndexInv(i),2+kk);
            if ConnList_global(nb,1) == 2
                num(1) = 78;%i;
            end
        end
    end
    if TipStates(IndexInv(i)) > 0.1 && AllEle_global(IndexInv(i),10)> 1.1 && ConnList(i,1) == 2
        num(2) = i;
    end
%     if TipStates(IndexInv(i)) > 0.1 
%         num(1) = i;
%     end
    if AllEle_global(IndexInv(i),10)> 1.1 && ConnList(i,1) == 2
        if DD(i+nEle) < -1e-15
            for kk = 1 : ConnList_global(IndexInv(i),2)
                nb = ConnList_global(IndexInv(i),2+kk);
                if Index(nb) > 0.1 && abs(DD_global(MaxEle+nb)) < 1e-20 && ConnList_global(nb,2) < 2.1
                    num(3) = i;
                end
            end
        end
        
    end
    
end
disp('Show Displacement and Stress around Fracture Surface');
disp('');
eps = 0.001;
% Draw NF
TangS = zeros(3,1);% Tangential for 3 tips
%% KI at the intersection
% Search
a = -120;
b = 120;
if num(1) > 0.1
    L = AllEle(num(1),7);
    p0 = [AllEle(num(1),8),AllEle(num(1),9)];
    for i = 1 : ConnList_global(IndexInv(num(1)),2)
        nb = ConnList_global(IndexInv(num(1)),2+i);
        if nb > 0.1 && AllEle_global(nb,10) < 1.1
            p1 = [AllEle_global(nb,8),AllEle_global(nb,9)];
        end
    end
    dis = p0-p1;
    vec = dis/sqrt(dis(1)^2+dis(2)^2);
    vec0 = [AllEle(num(1),6) AllEle(num(1),5)];
    if dot(vec,vec0) > 1e-8
        %正向
        type = 1;
    else
        type = 2;
    end
    Sn = -1E10;
    Tauv = 1E10;
    thetai= 0;
    sigmaN = zeros(200,1);
    tau = zeros(200,1);
    count = 0;
    for i = a:b
        count = count + 1;
        alpha = asind(AllEle(num(1),5));
        if type < 1.1
            ppx = p0(1) + 0*L/2*AllEle(num(1),6);
            ppy = p0(2) + 0* L/2*AllEle(num(1),5);
            xi = ppx + eps*cosd(i+alpha);
            yi = ppy + eps*sind(i+alpha);
            [~,~,Sx,Sy,Sxy,~,~] =  CalcPointStress_C_local(Mat,DD,xi,yi,AllEle);
            cosi = -sind(i+alpha);
            sini = cosd(i+alpha);
            Sxxi = Sx + Mat.Sxx;
            Syyi = Sy + Mat.Syy;
            Sxyi = Sxy + Mat.Sxy;
            sigmaN(count) = (Sxxi * cosi^2 + 2 * Sxyi * sini * cosi + Syyi * sini^2);
            tau(count) = (Sxxi- Syyi)*sini*cosi-Sxyi*(cosi^2 - sini^2);
            if sigmaN(count)> Sn
                Sn = sigmaN(count);
                Tauv = tau(count);
                thetai = i;
            end
        else
            ppx = p0(1) - 0* L/2*AllEle(num(1),6);
            ppy = p0(2) -  0*L/2*AllEle(num(1),5);
            xi = ppx + eps*cosd(-i+alpha-180);
            yi = ppy + eps*sind(-i+alpha-180);
            [~,~,Sx,Sy,Sxy,~,~] =  CalcPointStress_C_local(Mat,DD,xi,yi,AllEle);
            Sxxi = Sx+ Mat.Sxx;
            Syyi = Sy + Mat.Syy;
            Sxyi = Sxy + Mat.Sxy;
            cosi = -sind(-i+alpha-180);
            sini = cosd(-i+alpha-180);
            sigmaN(count) = (Sxxi * cosi^2 + 2 * Sxyi * sini * cosi + Syyi * sini^2);
            tau(count) = (Sxxi- Syyi)*sini*cosi-Sxyi*(cosi^2 - sini^2);
            if sigmaN(count) > Sn
                Sn = sigmaN(count);
                Tauv = tau(count);
                thetai = i;
            end
        end
    end
    K11 = Sn;%*sqrt(2*pi*eps);
    theta11 = thetai;
%% Find Critial length of tensil stess
To = 1;% 1 MPa
HasFind = 0;
tol = 0.1;
Lleft = 1e-5;
Lright = 1e-1;
eps = 0.5*(Lleft+Lright);
while HasFind < 0.1
    if type < 1.1
        ppx = p0(1) + L/2*AllEle(num(1),6);
        ppy = p0(2) + L/2*AllEle(num(1),5);
        xi = ppx + eps*cosd(theta11+alpha);
        yi = ppy + eps*sind(theta11+alpha);
        [~,~,Sx,Sy,Sxy,~,~] =  CalcPointStress_C_local(Mat,DD,xi,yi,AllEle);
        cosi = -sind(theta11+alpha);
        sini = cosd(theta11+alpha);
        Sxxi = Sx + Mat.Sxx;
        Syyi = Sy + Mat.Syy;
        Sxyi = Sxy + Mat.Sxy;
        sigmaN= (Sxxi * cosi^2 + 2 * Sxyi * sini * cosi + Syyi * sini^2);
    else
        ppx = p0(1) - L/2*AllEle(num(1),6);
        ppy = p0(2) - L/2*AllEle(num(1),5);
        xi = ppx + eps*cosd(-theta11+alpha-180);
        yi = ppy + eps*sind(-theta11+alpha-180);
        [~,~,Sx,Sy,Sxy,~,~] =  CalcPointStress_C_local(Mat,DD,xi,yi,AllEle);
        Sxxi = Sx+ Mat.Sxx;
        Syyi = Sy + Mat.Syy;
        Sxyi = Sxy + Mat.Sxy;
        cosi = -sind(-theta11+alpha-180);
        sini = cosd(-theta11+alpha-180);
        sigmaN = (Sxxi * cosi^2 + 2 * Sxyi * sini * cosi + Syyi * sini^2);
    end
    [Lleft,Lright,eps,HasFind] = FindCritDisT(Lleft,Lright,eps,sigmaN,To,tol);
end
% Calculated Energy
G = 0;
%  integration
Lc = eps;
xv = 1e-4:1e-4:Lc; 
for i = 1 : length(xv)
    eps = xv(i);
    if type < 1.1
        ppx = p0(1) + L/2*AllEle(num(1),6);
        ppy = p0(2) + L/2*AllEle(num(1),5);
        xi = ppx + eps*cosd(theta11+alpha);
        yi = ppy + eps*sind(theta11+alpha);
        [~,~,Sx,Sy,Sxy,~,~] =  CalcPointStress_C_local(Mat,DD,xi,yi,AllEle);
        cosi = -sind(theta11+alpha);
        sini = cosd(theta11+alpha);
        Sxxi = Sx + Mat.Sxx;
        Syyi = Sy + Mat.Syy;
        Sxyi = Sxy + Mat.Sxy;
        sigmaN= (Sxxi * cosi^2 + 2 * Sxyi * sini * cosi + Syyi * sini^2);
    else
        ppx = p0(1) - L/2*AllEle(num(1),6);
        ppy = p0(2) - L/2*AllEle(num(1),5);
        xi = ppx + eps*cosd(-theta11+alpha-180);
        yi = ppy + eps*sind(-theta11+alpha-180);
        [~,~,Sx,Sy,Sxy,~,~] =  CalcPointStress_C_local(Mat,DD,xi,yi,AllEle);
        Sxxi = Sx+ Mat.Sxx;
        Syyi = Sy + Mat.Syy;
        Sxyi = Sxy + Mat.Sxy;
        cosi = -sind(-theta11+alpha-180);
        sini = cosd(-theta11+alpha-180);
        sigmaN = (Sxxi * cosi^2 + 2 * Sxyi * sini * cosi + Syyi * sini^2);
    end
    G = G + (sigmaN*sqrt(eps*pi/2))^2;
end
KC = sqrt(G);
    
%% KI at boundary between OpenNF and closeNF
% Search
if num(3) > 0.1
    p0 = [AllEle(num(3),8),AllEle(num(3),9)];
    for i = 1 : ConnList_global(IndexInv(num(3)),2)
        nb = ConnList_global(IndexInv(num(3)),2+i);
        if nb > 0.1 && abs(DD_global(MaxEle+nb)) < 1e-20 && ConnList_global(nb,2) < 2.1
            p1 = [AllEle_global(nb,8),AllEle_global(nb,9)];
        end
    end
    dis = p0-p1;
    vec = dis/sqrt(dis(1)^2+dis(2)^2);
    vec0 = [AllEle(num(3),6) AllEle(num(3),5)];
    if dot(vec,vec0) > 1e-8
        %正向
        type = 2;
    else
        type = 1;
    end
    Sn = -1E10;
    Tauv = 1E10;
    thetai= 0;
    L = AllEle(num(3),7);
    sigmaN = zeros(200,1);
    tau = zeros(200,1);
    count = 0;
    for i = a:b
        count = count + 1;
        alpha = asind(AllEle(num(3),5));
        if type < 1.1
            ppx = p0(1) + L/2*AllEle(num(3),6);
            ppy = p0(2) + L/2*AllEle(num(3),5);
            xi = ppx + eps*cosd(i+alpha);
            yi = ppy + eps*sind(i+alpha);
            [~,~,Sx,Sy,Sxy,~,~] =  CalcPointStress_C_local(Mat,DD,xi,yi,AllEle);
            cosi = -sind(i+alpha);
            sini = cosd(i+alpha);
            Sxxi = Sx+ Mat.Sxx;
            Syyi = Sy + Mat.Syy;
            Sxyi = Sxy + Mat.Sxy;
            sigmaN(count) = (Sxxi * cosi^2 + 2 * Sxyi * sini * cosi + Syyi * sini^2);
            tau(count) = (Sxxi- Syyi)*sini*cosi-Sxyi*(cosi^2 - sini^2);
            if sigmaN(count)> Sn
                Sn = sigmaN(count);
                Tauv = tau(count);
                thetai = i;
            end
        else
            ppx = p0(1) - L/2*AllEle(num(3),6);
            ppy = p0(2) - L/2*AllEle(num(3),5);
            xi = ppx + eps*cosd(-i+alpha-180);
            yi = ppy + eps*sind(-i+alpha-180);
            [~,~,Sx,Sy,Sxy,~,~] =  CalcPointStress_C_local(Mat,DD,xi,yi,AllEle);
            Sxxi = Sx+ Mat.Sxx;
            Syyi = Sy + Mat.Syy;
            Sxyi = Sxy + Mat.Sxy;
            cosi = -sind(-i+alpha-180);
            sini = cosd(-i+alpha-180);
            sigmaN(count) = (Sxxi * cosi^2 + 2 * Sxyi * sini * cosi + Syyi * sini^2);
            tau(count) = (Sxxi- Syyi)*sini*cosi-Sxyi*(cosi^2 - sini^2);
            if sigmaN(count) > Sn
                Sn = sigmaN(count);
                Tauv = tau(count);
                thetai = i;
            end
        end
    end
    K22 = Sn*sqrt(2*pi*eps);
    theta22 = thetai;
end
%% KI at shearing Tip
% Search
if num(2) > 0.1
    p0 = [AllEle(num(2),8),AllEle(num(2),9)];
    for i = 1 : ConnList_global(IndexInv(num(2)),2)
        nb = ConnList_global(IndexInv(num(2)),2+i);
        if Index(nb) > 0.1
            p1 = [AllEle_global(nb,8),AllEle_global(nb,9)];
        end
    end
    dis = p0-p1;
    vec = dis/sqrt(dis(1)^2+dis(2)^2);
    vec0 = [AllEle(num(2),6) AllEle(num(2),5)];
    if dot(vec,vec0) > 1e-8
        %正向
        type = 1;
    else
        type = 2;
    end
    Sn = -1E10;
    Tauv = 1E10;
    thetai= 0;
    L = AllEle(num(2),7);
    sigmaN = zeros(200,1);
    tau = zeros(200,1);
    count = 0;
    for i = a:b
        count = count + 1;
        alpha = asind(AllEle(num(2),5));
        if type < 1.1
            ppx = p0(1) + L/2*AllEle(num(2),6);
            ppy = p0(2) + L/2*AllEle(num(2),5);
            xi = ppx + eps*cosd(i+alpha);
            yi = ppy + eps*sind(i+alpha);
            [~,~,Sx,Sy,Sxy,~,~] =  CalcPointStress_C_local(Mat,DD,xi,yi,AllEle);
            cosi = -sind(i+alpha);
            sini = cosd(i+alpha);
            Sxxi = Sx+ Mat.Sxx;
            Syyi = Sy + Mat.Syy;
            Sxyi = Sxy + Mat.Sxy;
            sigmaN(count) = (Sxxi * cosi^2 + 2 * Sxyi * sini * cosi + Syyi * sini^2);
            tau(count) = (Sxxi- Syyi)*sini*cosi-Sxyi*(cosi^2 - sini^2);
            if sigmaN(count)> Sn
                Sn = sigmaN(count);
                Tauv = tau(count);
                thetai = i;
            end
        else
            ppx = p0(1) - L/2*AllEle(num(2),6);
            ppy = p0(2) - L/2*AllEle(num(2),5);
            xi = ppx + eps*cosd(-i+alpha-180);
            yi = ppy + eps*sind(-i+alpha-180);
            [~,~,Sx,Sy,Sxy,~,~] =  CalcPointStress_C_local(Mat,DD,xi,yi,AllEle);
            Sxxi = Sx+ Mat.Sxx;
            Syyi = Sy + Mat.Syy;
            Sxyi = Sxy + Mat.Sxy;
            cosi = -sind(-i+alpha-180);
            sini = cosd(-i+alpha-180);
            sigmaN(count) = (Sxxi * cosi^2 + 2 * Sxyi * sini * cosi + Syyi * sini^2);
            tau(count) = (Sxxi- Syyi)*sini*cosi-Sxyi*(cosi^2 - sini^2);
            if sigmaN(count) > Sn
                Sn = sigmaN(count);
                Tauv = tau(count);
                thetai = i;
            end
        end
    end
    K33 = Sn*sqrt(2*pi*eps);
    theta33 = thetai;
end
clear DD AllEle ConnList

end

