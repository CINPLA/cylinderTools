function f = plot_points_and_cylinder( points, cylinder )
%PLOT_POINTS_ Summary of this function goes here
%   Detailed explanation goes here

M = cylinder.M;
r0 = cylinder.r0; 
p = cylinder.p;

points(1,:) = points(1,:) - r0(1);
points(2,:) = points(2,:) - r0(2);
points(3,:) = points(3,:) - r0(3);
points = M*points;

xp = points(1,:); 
yp = points(2,:); 
zp = points(3,:);

f = figure;
scatter3(xp, yp, zp, 'filled')



N = 100;
d = 20;
xx = linspace(min(xp)-d, max(xp)+d, N);
yy = linspace(min(yp)-d, max(yp)+d, N);
zz = linspace(min(zp)-d, max(zp)+d, N); 

[xx,yy,zz] = meshgrid(xx,yy,zz); 

F = (p(1)*xx.*xx + p(3)*yy.*yy + p(2)*xx.*yy + p(4)*xx + p(5)*yy + p(6));
isosurface(xx,yy,zz,F,0);
axis equal
title('3d reconstruction of muscle fiber with nuclei locations')
xlabel('x position [um]')
ylabel('y position [um]')
zlabel('z position [um]')

end



