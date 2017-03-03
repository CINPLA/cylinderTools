function cylinder = fit_cylinder( points )
%FIT_CYLINDER Summary of this function goes here
%   Detailed explanation goes here

r0 = [mean(points(1,:)), mean(points(2,:)), mean(points(3,:))];

points(1,:) = points(1,:) - r0(1);
points(2,:) = points(2,:) - r0(2);
points(3,:) = points(3,:) - r0(3);

lb = [-pi, -pi, -pi/2];
ub = [pi, pi, pi/2]; 

x0 = [0 0 0];
myfun = @(x) cylinderDistance(x, points);

function [c,ceq] = mycon(x)
    Mcon = ang2rot(x(1), x(2), x(3)); 
    pointscon = Mcon*points;
    xcon = pointscon(1,:); 
    ycon = pointscon(2,:);
    p = ellipsefit_direct(xcon,ycon);
    a = p(1); 
    b = p(2)/2;
    c = p(3); 
    d = p(4)/2; 
    f = p(5)/2;
    g = p(6);
    numerator = 2*(a*f*f + c*d*d + g*b*b - 2*b*d*f - a*c*g); 
    denom_sqrt = sqrt((a - c)^2 + 4*b*b);
    denominator_a = (b*b - a*c)*(denom_sqrt - (a+c));
    denominator_b = (b*b - a*c)*(-denom_sqrt - (a+c));

	a = sqrt(numerator/denominator_a);
	b = sqrt(numerator/denominator_b); 
    
    aa = max([a,b]);
    bb = min([a,b]);
    
    c = aa-2.0*bb;
    ceq = [];
end

x = fmincon(myfun,x0,[],[],[],[],lb,ub, @mycon);
M = ang2rot(x(1), x(2), x(3)); 
% 
points_aligned = M*points;
% 
% 
% 
xp = points_aligned(1,:); 
yp = points_aligned(2,:);
% zp = points_aligned(3,:);
% 
% 
% figure
% hold on
% scatter3(xp, yp, zp, 'r')
% 
p = ellipsefit_direct(xp,yp);

cylinder = {}; 

cylinder.M = M;
cylinder.angles = [x(1), x(2), x(3)];
cylinder.r0 = r0;
cylinder.p = p;


a = p(1); 
b = p(2)/2;
c = p(3); 
d = p(4)/2; 
f = p(5)/2;
g = p(6);

cylinder.ellipse = {};
cylinder.ellipse.x0 = (c*d - b*f)/(b*b - a*c);
cylinder.ellipse.y0 = (a*f - b*d)/(b*b - a*c); 

numerator = 2*(a*f*f + c*d*d + g*b*b - 2*b*d*f - a*c*g); 
denom_sqrt = sqrt((a - c)^2 + 4*b*b);
denominator_a = (b*b - a*c)*(denom_sqrt - (a+c));
denominator_b = (b*b - a*c)*(-denom_sqrt - (a+c));

cylinder.ellipse.a = sqrt(numerator/denominator_a);
cylinder.ellipse.b = sqrt(numerator/denominator_b); 


if b == 0
    if a < c
        phi = 0;
    else
        phi = pi/2; 
    end
else
    if a < c
        phi = 0.5*acot((a-c)/(2*b));
    else
        phi = pi/2 + 0.5*acot((a-c)/(2*b));
    end
end

cylinder.ellipse.phi = phi;

if cylinder.ellipse.a < cylinder.ellipse.b
    a = cylinder.ellipse.a;
    b = cylinder.ellipse.b;
    cylinder.ellipse.a = b; 
    cylinder.ellipse.b = a;
    cylinder.ellipse.phi = cylinder.ellipse.phi + pi/2;
end


%cylinder.ellipse.a = sqrt( / ((b*b - a*c)*(sqrt(a))))






% 
% N = 100;
% d = 20;
% xx = linspace(min(xp)-d, max(xp)+d, N);
% yy = linspace(min(yp)-d, max(yp)+d, N);
% zz = linspace(min(zp)-d, max(zp)+d, N); 
% 
% [xx,yy,zz] = meshgrid(xx,yy,zz); 
% F = (p(1)*xx.*xx + p(3)*yy.*yy + p(2)*xx.*yy + p(4)*xx + p(5)*yy + p(6));
% isosurface(xx,yy,zz,F,0);
% axis equal
end

