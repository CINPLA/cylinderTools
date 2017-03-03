function avgnnd = random_avg_nnd(cylinder_projections)
%RANDOM_AVG_NND Summary of this function goes here
%   Detailed explanation goes here

nExperiments = 10000;

L = cylinder_projections.zmax - cylinder_projections.zmin;
N = length(cylinder_projections.z_position);

nd_all = [];

for i=1:nExperiments
    cnew = cylinder_projections;
    zpos = rand(1,N)*L + cylinder_projections.zmin;
    cnew.z_position = zpos; 
    arcpos = rand(1,N)*cylinder_projections.arclength;
    cnew.arc_position = arcpos; 
    
    [ neighbor_idx, nearest_dist] = periodic_nearest_neighbor( cnew );
    nd_all = [nd_all, nearest_dist];
    
end

avgnnd = mean(nd_all);

end

