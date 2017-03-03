function radial_distribution = radial_distribution( cylinder_projections )
%RADIAL_DISTRIBUTION Summary of this function goes here
%   Detailed explanation goes here
zp = cylinder_projections.z_position;
zlen = cylinder_projections.zmax - cylinder_projections.zmin;
ap = cylinder_projections.arc_position;
alen = cylinder_projections.arclength;
z_pos = [zp-zlen, zp-zlen, zp-zlen, zp, zp, zp, zp+zlen, zp+zlen, zp+zlen];
a_pos = [ap-alen, ap, ap+alen, ap-alen, ap, ap+alen, ap-alen, ap, ap+alen];





pos_periodic = [z_pos; a_pos];
pos = [zp; ap];
pd = pdist2(pos', pos_periodic');
%f = figure;
pd = pd(:); 
pd(pd==0) = [];
pd(pd>200) = [];
h = histogram(pd);
radial_distribution = {}; 
radial_distribution.pairwise_distance = pd; 
radial_distribution.radial_distribution_histogram = h;
%radial_distribution.figure = f;
end

