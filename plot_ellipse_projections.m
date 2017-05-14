function f = plot_ellipse_projections( points, cylinder )
%PLOT_ELLIPSE_PROJECTIONS Summary of this function goes here
%   Detailed explanation goes here

f = figure;
M = cylinder.M;
r0 = cylinder.r0; 
p = cylinder.p;
ellipse = cylinder.ellipse;


points(1,:) = points(1,:) - r0(1);
points(2,:) = points(2,:) - r0(2);
points(3,:) = points(3,:) - r0(3);
points = M*points;


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

[~,~,xf,yf] = ellipse_distance(xp, yp, p);

% plot_points_and_ellipse(points, cylinder);



for i=1:length(xp)
    plot([xp(i), xf(i)], [yp(i), yf(i)], 'b-');
end
axis equal

title('Nuclei projected on fiber surface')
xlabel('x position [um]')
ylabel('y position [um]')


end
