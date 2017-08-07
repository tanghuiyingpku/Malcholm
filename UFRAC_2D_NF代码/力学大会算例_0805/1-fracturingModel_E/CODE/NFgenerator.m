function NFgenerator()
theta0 = -15;
% Horizontal wells in the middle
Length = 0.1;
SpaceX = 0.1;
SpaceY = 0.08;
Lx = 1;
Ly = 0.8;

xsc = 50;
ysc = 50;

SpaceX = SpaceX*xsc;
SpaceY = SpaceY*ysc;
Lx = Lx*xsc;
Ly = Ly*ysc;
Length = Length * xsc;

nX = floor(Lx/SpaceX);
nY = floor(Ly/SpaceY);
y0 =0;%1000;
NFcenter = zeros(nX*nY,2); %[x,y]
nF = 0;
iy = 0;
while iy < nY-1
    iy = iy + 1;
    if mod(iy,2) < 0.1
        nF0 = nF;
        nF = nF + nX-2;
        y = iy * SpaceY+y0;
        x = SpaceX*1.5:SpaceX:Lx-SpaceX*1.5;
        %x = x + randn()*50;
        NFcenter(nF0+1:nF,:) = [x',ones(length(x),1)*y];
    else
        nF0 = nF;
        nF = nF + nX-1;
        y = iy * SpaceY+y0;
        x = SpaceX:SpaceX:Lx-SpaceX;
        %x = x + randn()*50;
        NFcenter(nF0+1:nF,:) = [x',ones(length(x),1)*y];
    end
end
NFnet = zeros(nF,4);% x1 y1 x2 y2
figure(1);
for i =1 : nF
    theta = theta0 + randn()*0;
    NFnet(i,:) = [NFcenter(i,1)-Length/2*cosd(theta),NFcenter(i,2)-Length/2*sind(theta),NFcenter(i,1)+Length/2*cosd(theta),NFcenter(i,2)+Length/2*sind(theta)];
    plot([NFnet(i,1),NFnet(i,3)],[NFnet(i,2),NFnet(i,4)]);
    hold on;
end
% NFnet(:,1)= NFnet(:,1)/Lx;
% NFnet(:,3)= NFnet(:,3)/Lx;
% NFnet(:,2)= NFnet(:,2)/Ly;
% NFnet(:,4)= NFnet(:,4)/Ly;
%axis([0 2500 0 2000]);
axis equal;
save GenerateFrac.frac -ascii NFnet;
% %another set
% theta0 = -30;
% % Horizontal wells in the middle
% Length = 180;
% SpaceX = 200;
% SpaceY = 250;
% Lx = 750;
% Ly = 1500;
% nX = floor(Lx/SpaceX);
% nY = floor(Ly/SpaceY);
% y0 =0;%1000;
% NFcenter = zeros(nX*nY,2); %[x,y]
% nF = 0;
% iy = 0;
% while iy < nY-1
%     iy = iy + 1;
%     if mod(iy,2) > 0.1
%         nF0 = nF;
%         nF = nF + nX-2;
%         y = (iy + 0.5) * SpaceY+y0;
%         x = SpaceX*1.5:SpaceX:Lx-SpaceX*1.5;
%         %x = x + randn()*50;
%         NFcenter(nF0+1:nF,:) = [x',ones(length(x),1)*y];
%     else
%         nF0 = nF;
%         nF = nF + nX-1;
%         y = (iy + 0.5)* SpaceY+y0;
%         x = SpaceX:SpaceX:Lx-SpaceX;
%         %x = x + randn()*50;
%         NFcenter(nF0+1:nF,:) = [x',ones(length(x),1)*y];
%     end
% end
% NFnet = zeros(nF,5);% x1 y1 x2 y2
% figure(1);
% for i =1 : nF
%     theta = theta0 + randn()*25;
%     NFnet(i,:) = [NFcenter(i,1)-Length/2*cosd(theta),NFcenter(i,2)-Length/2*sind(theta),NFcenter(i,1)+Length/2*cosd(theta),NFcenter(i,2)+Length/2*sind(theta),0];
%     plot([NFnet(i,1),NFnet(i,3)],[NFnet(i,2),NFnet(i,4)]);
%     hold on;
% end
% NFnet(:,1)= NFnet(:,1)/Lx;
% NFnet(:,3)= NFnet(:,3)/Lx;
% NFnet(:,2)= NFnet(:,2)/Ly;
% NFnet(:,4)= NFnet(:,4)/Ly;
% axis([0 2500 0 2000]);
% axis equal;
% save GenerateFrac2.frac -ascii NFnet;
f = 1;




