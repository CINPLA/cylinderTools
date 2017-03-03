function arc_length = NewArcLengthEllipse( ellipse, phi_1, phi_0 )
%NEWARCLENGTHELLIPSE Summary of this function goes here
%   Detailed explanation goes here

a = ellipse.a; 
b = ellipse.b; 
s = @(x) abs(sqrt(1+(x.*x*b*b/(a*a))./(a*a-x.*x)));

s_half = integral(s, -a, a);

s_list = [];
for phi=[phi_0, phi_1]
    n = floor(phi/pi);
    arc_length = n*s_half;
        
    phi = phi-pi*n;
    s_list(end+1) = n*abs(s_half) + abs(integral(s,min(a*cos(phi), a-1e-8), a-1e-8));
end
    
arc_length = s_list(2)-s_list(1);   
end

