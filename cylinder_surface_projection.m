function cylinder_projection = cylinder_surface_projection( points, cylinder )
%CYLINDER_SURFACE_PROJECTION Summary of this function goes here
%   Detailed explanation goes here
M = cylinder.M;
r0 = cylinder.r0; 
p = cylinder.p;
ellipse = cylinder.ellipse;


points(1,:) = points(1,:) - r0(1);
points(2,:) = points(2,:) - r0(2);
points(3,:) = points(3,:) - r0(3);
points = M*points;

points2d = [points(1,:) - ellipse.x0; points(2,:) - ellipse.y0];


R = [cos(ellipse.phi), -sin(ellipse.phi); sin(ellipse.phi), cos(ellipse.phi)];

points2d = inv(R)*points2d;

z_position = points(3,:);
x_pos = points(1,:);
y_pos = points(2,:);
[~,~,xf,yf] = ellipse_distance(x_pos, y_pos, p);
scatter(x_pos, y_pos)

arc_position = zeros(size(z_position));
theta_list = zeros(size(z_position));
arclength = NewArcLengthEllipse(ellipse, 2*pi,0);
for i=1:length(points(1,:))
    theta = atan2(yf(i),xf(i));
    theta_list(i) = theta;
    arc_position(i) = NewArcLengthEllipse(ellipse, theta, 0);
    if arc_position(i) < 0
        arc_position(i) = arc_position(i) + arclength;
    end
end

d = 10;
zmin = min(z_position) - d;
zmax = max(z_position) + d;
arclength = NewArcLengthEllipse(ellipse, 2*pi,0);

cylinder_projection = {};
cylinder_projection.zmin = zmin;
cylinder_projection.zmax = zmax;
cylinder_projection.arclength = arclength;
cylinder_projection.z_position = z_position;
cylinder_projection.arc_position = arc_position;
cylinder_projection.x_position = xf;
cylinder_projection.y_position = yf;
cylinder_projection.theta = theta_list;

end

