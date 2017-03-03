function voro3d = voronoi3d( cylinder_struct, points )
%VORONOI3D Summary of this function goes here
%   Detailed explanation goes here

M = cylinder_struct.M;
r0 = cylinder_struct.r0; 
p = cylinder_struct.p;
ellipse = cylinder_struct.ellipse;

points(1,:) = points(1,:) - r0(1);
points(2,:) = points(2,:) - r0(2);
points(3,:) = points(3,:) - r0(3);
points = M*points;


points2d = [points(1,:) - ellipse.x0; points(2,:) - ellipse.y0];
R = [cos(ellipse.phi), -sin(ellipse.phi); sin(ellipse.phi), cos(ellipse.phi)];
points2d = inv(R)*points2d;
points(1:2,:) = points2d;

philist = linspace(0,2*pi,100); 
xlist = ellipse.a*cos(philist); 
ylist = ellipse.b*sin(philist);
xvec = [xlist; ylist];

zmin = min(points(3,:));
zmax = max(points(3,:));
d = 0.1*(zmax - zmin);
zmin = zmin-d; 
zmax = zmax+d;

nPoints3d = 200;
x = linspace(-ellipse.a*1.1, ellipse.a*1.1,nPoints3d);
y = linspace(-ellipse.b*1.1, ellipse.b*1.1,nPoints3d);
z = linspace(zmin,zmax,nPoints3d);
[X,Y,Z] = meshgrid(x,y,z);
xmin = min(X(:));
ymin = min(Y(:));
xmax = max(X(:));
ymax = max(Y(:));

[Xe,Ye,Ze] = cylinder(linspace(1,1,nPoints3d),nPoints3d);
Ze = Ze*(zmax-zmin) + zmin;
Xe = Xe*ellipse.a;
Ye = Ye*ellipse.b; 

[Xt,Yt,Zt] = sphere(nPoints3d);
Zt(:,:) = zmax;
Xt = Xt*ellipse.a;
Yt = Yt*ellipse.b; 
 
[Xb,Yb,Zb] = sphere(nPoints3d);
Zb(:,:) = zmin;
Xb = Xb*ellipse.a;
Yb = Yb*ellipse.b; 

[IDX3D, dist3D] = knnsearch(points', [X(:), Y(:), Z(:)]);
IDX3D = reshape(IDX3D, nPoints3d, nPoints3d, nPoints3d);
dist3D = reshape(dist3D, nPoints3d, nPoints3d, nPoints3d);
newidx3d = randperm(max(IDX3D(:)));
IDX3Dn = zeros(size(IDX3D)); 
for i=1:max(IDX3D(:))
    IDX3Dn(IDX3D==i) = newidx3d(i);
end
IDX3D = IDX3Dn;

volume_per_pixel = (max(X(:) - min(X(:))))*(max(Y(:) - min(Y(:))))*(max(Z(:) - min(Z(:))))/nPoints3d^3;
volumes = zeros(1, length(points(1,:)));
for i=1:length(points(1,:))
    volumes(i) = sum(IDX3D(:) == i)*volume_per_pixel;
end
f3 = figure;
boxplot(volumes)
title('Box plot of volumes')
ylabel('Volume (um^3)') 

f1=figure;
subplot(1,2,1)

slice(X,Y,Z,IDX3D, Xe, Ye, Ze);

hold on


slice(X,Y, Z, IDX3D,Xt, Yt, Zt )

slice(X,Y, Z, IDX3D,Xb, Yb, Zb )
axis equal
a = get(gca, 'children');
for i=1:length(a)
   set(a(i), 'edgecolor', 'none')
end

zlist = linspace(zmin, zmax, 10);
for i=1:length(zlist)
   plot3(xvec(1,:), xvec(2,:), linspace(zlist(i), zlist(i), 100), 'k')
end
title('3D Voronoi diagram')
xlabel('x [um]')
ylabel('y [um]')
zlabel('z [um]')
camlight headlight
lighting phong

subplot(1,2,2)
slice(X,Y,Z,dist3D, Xe, Ye, Ze);

hold on

slice(X,Y, Z, dist3D,Xt, Yt, Zt )

slice(X,Y, Z, dist3D,Xb, Yb, Zb )
axis equal
a = get(gca, 'children');
for i=1:length(a)
   set(a(i), 'edgecolor', 'none')
end
for i=1:length(zlist)
   plot3(xvec(1,:), xvec(2,:), linspace(zlist(i), zlist(i), 100), 'k')
end

camlight headlight
lighting phong
xlabel('x [um]')
ylabel('y [um]')
zlabel('z [um]')
title('3D distance map')

h = uicontrol('style','slider','units','pixel','position',[0 0 200 20]);
addlistener(h,'ContinuousValueChange',@(hObject, event) plotline_hObject_thing(hObject, event, ellipse,X, Y,Z, IDX3D,dist3D, f, 1:nPoints3d));


% nPoints3d = 200;
% zmin = min(points(3,:));
% zmax = max(points(3,:));
% d = 0.1*(zmax - zmin);
% zmin = zmin-d; 
% zmax = zmax+d;
% [Xe,Ye,Ze] = cylinder(linspace(1,1,nPoints3d),300);
% Ze = Ze*(zmax-zmin)+zmin;
% Xe = Xe*ellipse.a;
% Ye = Ye*ellipse.b; 
% 
% cs = cos(ellipse.phi); 
% sc = sin(ellipse.phi); 
% R = [cs, -sc; sc, cs];
% 
% points2d = points(1:2, :); 
% points2d = inv(R)*points2d;
% 
% points(1:2,:) = points2d;
% 
% scatter3(Xe(:), Ye(:), Ze(:))
% hold on
% scatter3(points(1,:), points(2,:), points(3,:))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f2 = figure;
% hold on
for i=1:nPoints3d
    for j=1:nPoints3d
        xvec = [x(i), y(j)];
        xvec = xvec';
        r = (xvec(1)/ellipse.a)^2 + (xvec(2)/ellipse.b)^2; 
        if r > 1
            IDX3D(i,j,:) = 0;
        end
    end
end

v_list = {}; 
f_list = {};
% Ftot = 
for i=1:max(IDX3D(:))
    IDX3Dtemp = IDX3D == i;
    IDX3Dtemp(:,:,1) = 0;
    IDX3Dtemp(:,:,end) = 0;
    if sum(IDX3Dtemp(:)) ~= 0
        [F,V] = MarchingCubes(X,Y,Z, IDX3Dtemp,0.3);
%        p = points(:,i);
%        p2 = p';
        col = [rand rand rand];
        patch('faces', F, 'vertices', V, 'facecolor', col, 'edgecolor', 'none')
%         hold on
    end
end
camlight right
lighting phong

voro3d = {};
voro3d.volumes = volumes; 
voro3d.slices = f1;
voro3d.volume_box_plot = f3; 
voro3d.fancy = f2;




% sbw = size(bw);
% [XX,YY,ZZ] = meshgrid(linspace(1,nImage,sbw(1)), linspace(1,mImage,sbw(2)), umPerSlice*(1:sbw(3)));
% [F,V] = MarchingCubes(XX,YY,ZZ, bw,0.3);
% 
% % for i=1:size(Xe,1)
% %     for j=1:size(Xe,2)
% %         xvec = [Xe(i,j)+ellipse.X0, Ye(i,j)+ellipse.Y0, Ze(i,j)];
% %         xvec = M*xvec';
% %         Xe(i,j) = xvec(1); 
% %         Ye(i,j) = xvec(2); 
% %         Ze(i,j) = xvec(3);
% %     end
% % end
% 
% [Xe,Ye,Ze] = cylinder(linspace(1,1,nPoints3d),300);
% Ze = Ze*500;
% Xe = Xe*ellipse.a;
% Ye = Ye*ellipse.b; 
% for i=1:length(Xe(:,1))
%     for j=1:length(Xe(1,:))
%         xvect = [Xe(i,j); Ye(i,j)]; 
%         xvect = S*xvect;
%         Xe(i,j) = xvect(1);
%         Ye(i,j) = xvect(2);
%     end
% end
% 
% V = (M'*V')';
% V(:,1) = V(:,1) - ellipse.X0;
% V(:,2) = V(:,2) - ellipse.Y0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% figure 
% surf(Xe,Ye,Ze)
% a = get(gca, 'children');
% for i=1:length(a)
%    set(a(i), 'edgecolor', 'none', 'facecolor', 'blue')
% end
% alpha(0.5)
% hold on
% patch('faces', F, 'vertices', V, 'facecolor', 'r', 'edgecolor', 'none')
% camlight right
% lighting phong
% xlabel('x [um]')
% ylabel('y [um]')
% zlabel('z [um]')
% title ('3D reconstruction of cylinder')

end

