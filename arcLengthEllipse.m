function L = arcLengthEllipse( ellipse, phi_1, phi_0 )
%ARCLENGTHELLIPSE Summary of this function goes here
%   Detailed explanation goes here
a = ellipse.a; 
b = ellipse.b; 
e = sqrt(1-b^2/a^2);

t_1 = atan(a*tan(phi_1)/b);
t_0 = atan(a*tan(phi_0)/b);


L = b*(ellipticE(t_1,e) - ellipticE(t_0,e));

end

