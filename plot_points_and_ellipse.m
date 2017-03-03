function f = plot_points_and_ellipse( points, cylinder )
%PLOT_ELLIPSE Summary of this function goes here
%   Detailed explanation goes here

f = figure;

r0 = cylinder.r0;
points(1,:) = points(1,:) - r0(1);
points(2,:) = points(2,:) - r0(2);
points(3,:) = points(3,:) - r0(3);

points = cylinder.M*points;
ellipse = cylinder.ellipse; 


xp = points(1,:); 
yp = points(2,:);

scatter(xp, yp)

hold on

rho = linspace(0,2*pi, 1001); 

x = ellipse.a*cos(rho); 
y = ellipse.b*sin(rho); 

phi = ellipse.phi;
M = [cos(phi), -sin(phi); sin(phi), cos(phi)];

r = [x;y];
r = M*r;
xr = r(1,:); 
yr = r(2,:);

xr = xr + ellipse.x0;
yr = yr + ellipse.y0;

plot(xr,yr)
axis equal

end

