function [avg_3ddist, std_3ddist, avg_2ddist, std_2ddist, nnd] = optimal_distmaps(  cylinder, cylinder_projections, fname )
%OPTIMAL_DISTMAPS Summary of this function goes here
%   Detailed explanation goes here
n = length(cylinder_projections.x_position);
cylinder.M = eye(3); 
cylinder.angles = [0 0 0];
cylinder.r0 = [0 0 0];
cylinder.ellipse.x0 = 0; 
cylinder.ellipse.y0 = 0;
cylinder.ellipse.phi = 0;

Lx = cylinder_projections.zmax-cylinder_projections.zmin;
Ly = cylinder_projections.arclength;

xpos = cylinder_projections.z_position; 
ypos = cylinder_projections.arc_position;

xc = xpos; yc = ypos;

xpos = [xpos, xpos, xpos];
ypos = [ypos, ypos-Ly, ypos+Ly];

NN = 10000;



Nx = round(sqrt(Lx*NN/Ly));
Ny = round(NN/Nx);
NN = Nx*Ny;
x = linspace(cylinder_projections.zmin,cylinder_projections.zmax,Nx); 
y = linspace(0,Ly,Ny);
y = linspace(-Ly,2*Ly,3*Ny); 
[X,Y] = meshgrid(x,y);
grid = [X(:), Y(:)];



nIterations = 2000;
sigma = 10;

vx = zeros(1,n);
vy = zeros(1,n);

for i=1:nIterations

    points = [xpos; ypos]';
    
    [IDX, dist] = knnsearch(points,grid);
    IDX = reshape(IDX, 3*Ny,Nx);



    
    %IDX = mod(IDX, N);

    stats = regionprops(IDX+1,'Centroid'); 
    
    a = [stats.Centroid];
    a(isnan(a)) = [];
    a = reshape(a, 2,3*n);
    
    
    xc = a(1,1:n)/Nx*Lx + cylinder_projections.zmin;
    yc = mod(a(2,1:n),Ny)/Ny*Ly;
    
%     ax = n_xc - xc; 
%     ay = n_yc - yc;
%     
%     vx = 0.5*vx + ax; 
%     vy = 0.5*vy + ay;
%     
%     xc = xc + vx;
%     yc = yc + vy;
%     
%     xc = mod(xc, Lx)+cylinder_projections.zmin;
%     yc = mod(yc, Ly);
    
    xc = mod(xc + sigma*randn(1,n),Lx)+cylinder_projections.zmin;
    yc = mod(yc + sigma*randn(1,n),Ly);
    
    xpos = [xc, xc, xc];
    ypos = [yc, yc-Ly, yc+Ly];
    
    sigma = 0.99*sigma
%     scatter(xpos, ypos)
end

% x = linspace(cylinder_projections.zmin,cylinder_projections.zmax,Nx); 
% y = linspace(0, Ly,Ny); 
% [X,Y] = meshgrid(x,y);
% grid = [X(:), Y(:)];
% 
% [IDX, dist] = knnsearch(points,grid);
% IDX = reshape(IDX, Ny, Nx);
% IDX = mod(IDX, n);
% 
% imagesc(x,y,IDX)
% set(gca, 'ydir', 'normal')
% hold on
% scatter(xc, yc, 'filled', 'w')
% axis equal

new_cylinder_projections = cylinder_projections;
new_cylinder_projections.z_position = xc;
new_cylinder_projections.arc_position = yc;


t_len = 1000;
theta_table = linspace(0,2*pi-1e-9,t_len);
arclength_table = zeros(1,t_len); 
for i=1:t_len
    arclength_table(i) = NewArcLengthEllipse(cylinder.ellipse, theta_table(i), 0);
end


points = zeros(3,n);
for j=1:n
    points(3,j) = new_cylinder_projections.z_position(j);
    [m, idx] = min(abs(new_cylinder_projections.arc_position(j) - arclength_table));
    theta = theta_table(idx); 
    points(1,j) = cylinder.ellipse.a*cos(theta);
    points(2,j) = cylinder.ellipse.b*sin(theta);

end
new_cylinder_projections.x_position = points(1,:); 
new_cylinder_projections.y_position = points(2,:); 



voro2d = voronoi2d(new_cylinder_projections);

voro3d = voronoi3d( cylinder, points );

dd = voro3d.dist3d(:);
dd(isnan(dd)) = [];

avg_3ddist = mean(dd);
std_3ddist = std(dd); 
avg_2ddist = mean(voro2d.all_dist(:)); 
std_2ddist = std(voro2d.all_dist(:));
[ neighbor_idx, nnd] = periodic_nearest_neighbor(new_cylinder_projections );


figure(voro2d.voronoi)
print(voro2d.voronoi, '-dpdf', [fname '_voronoi2d.pdf']);

figure(voro2d.dist)
pdfname = [fname '_distance_map_2d.png'];
export_fig(pdfname, '-transparent');

figure(voro3d.slices)
pdfname = [fname '_distance_map_3d.png'];
export_fig(pdfname, '-transparent');

end

