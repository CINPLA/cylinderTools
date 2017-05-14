function [std_area, std_volumes, avg_3ddist, std_3ddist, avg_2ddist, std_2ddist] = random_distmaps( cylinder, cylinder_projections, nExperiments, fname)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n = length(cylinder_projections.x_position);

cylinder.M = eye(3); 
cylinder.angles = [0 0 0];
cylinder.r0 = [0 0 0];
cylinder.ellipse.x0 = 0; 
cylinder.ellipse.y0 = 0;
cylinder.ellipse.phi = 0;
arealist = [];
distlist2d = [];
vardistlist = [];
volumelist = [];
distlist3d = [];

t_len = 1000;
theta_table = linspace(0,2*pi-1e-9,t_len);
arclength_table = zeros(1,t_len); 
for i=1:t_len
    arclength_table(i) = NewArcLengthEllipse(cylinder.ellipse, theta_table(i), 0);
end




for i=1:nExperiments
    i
    random_cylinder_projections = cylinder_projections;
    random_cylinder_projections.arc_position = rand(1,n)*cylinder_projections.arclength;
    random_cylinder_projections.z_position = cylinder_projections.zmin + rand(1,n)*(cylinder_projections.zmax - cylinder_projections.zmin);
    
    points = zeros(3,n);
    for j=1:n
        points(3,j) = random_cylinder_projections.z_position(j);
        [m, idx] = min(abs(random_cylinder_projections.arc_position(j) - arclength_table));
        theta = theta_table(idx); 
        points(1,j) = cylinder.ellipse.a*cos(theta);
        points(2,j) = cylinder.ellipse.b*sin(theta);
        
    end
    voro2d = voronoi2d(random_cylinder_projections);
    arealist = [arealist, voro2d.arealist];
    distlist2d = [distlist2d, voro2d.all_dist(:)];
    vardistlist = [vardistlist, var(voro2d.all_dist)];
    %figure
    
    if i~=nExperiments
        close all
    end
    voro3d = voronoi3d( cylinder, points );
    volumelist = [volumelist, voro3d.volumes];
    dd = voro3d.dist3d(:);
    dd(isnan(dd)) = [];
    distlist3d = [distlist3d, dd];
    if i~=nExperiments
        close all
    end
    %hist(voronoi.all_dist)

    
end
std_area = std(arealist(:));
std_volumes = std(volumelist(:)); 
avg_3ddist = mean(distlist3d(:)); 
std_3ddist = std(distlist3d(:)); 
avg_2ddist = mean(distlist2d(:)); 
std_2ddist = std(distlist2d(:));



figure(voro2d.voronoi)
print(voro2d.voronoi, '-dpdf', [fname '_voronoi2d.pdf']);

figure(voro2d.dist)
pdfname = [fname '_distance_map_2d.png'];
export_fig(pdfname, '-transparent');

figure(voro3d.slices)
pdfname = [fname '_distance_map_3d.png'];
export_fig(pdfname, '-transparent');

% histogram(arealist)


end

