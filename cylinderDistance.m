function d = cylinderDistance( x, points )
%CYLINDERDISTANCE Summary of this function goes here
%   Detailed explanation goes here

M = ang2rot(x(1), x(2), x(3));

dvec = M*points;
dvec2 = dvec(1:2,:);
x = dvec(1,:);
y = dvec(2,:);
%ellip = fit_ellipse(dvec(:,1), dvec(:,2));

p = ellipsefit_direct(x,y);

[e,d,xf,yf] = ellipse_distance(x, y, p);

%ParG = [ellip.X0, ellip.Y0, ellip.a, ellip.b, ellip.phi];
%[sd, XYproj] = Residuals_ellipse(dvec2,ParG);

d = norm(d);
%d = sd

end

