function MapSurfS (S)
% This function is generating triangle maps
fid = fopen( ['..\TriangleGridding\3Tri2Grid\', 'dfm2D.grid' ] ,'r');
pointnum = fscanf(fid,'%d',1);
edgenum = fscanf(fid,'%d',1);
elenum =fscanf(fid,'%d',1);
temp = fgetl(fid);
point = zeros(pointnum,2);
for i = 1 : pointnum
    point(i,1:2) = fscanf(fid,'%f',2);
    temp = fgetl(fid);
end
fprintf('%s',temp);
edge = zeros(edgenum,2);
for i = 1 : edgenum
    edge(i,1:2) = fscanf(fid,'%f',2);
    temp = fgetl(fid);
end
fprintf('%s',temp);
ele = zeros(elenum,3);
Center = zeros(elenum,2);

for i = 1 : elenum
    temp = fscanf(fid,'%f',1);
    ele(i,1:3) = fscanf(fid,'%f',3);
    ele(i,1:3) = ele(i,1:3) + 1;
    temp = fgetl(fid);
end
fclose(fid);
% CalcCenter
global PicScale;
figure 
Pm = S;
for i = 1 : elenum
    px = [ point(ele(i,1),1), point(ele(i,2),1),point(ele(i,3),1) ];
    py = [ point(ele(i,1),2), point(ele(i,2),2),point(ele(i,3),2) ];
    Center(i,1) = sum(px)/3 - (PicScale(2) - PicScale(1))/2;
    Center(i,2) = sum(py)/3- (PicScale(4) - PicScale(3))/2;
end
[XX,YY] = meshgrid(linspace(PicScale(1),PicScale(2),50),linspace(PicScale(3),PicScale(4),50));
Pm2 = griddata(Center(:,1),Center(:,2),Pm,XX,YY);
surf(XX,YY,Pm2);
shading flat;
view(2);
end
