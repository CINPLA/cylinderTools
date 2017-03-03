function  [ neighbor_idx, nearest_dist] = periodic_nearest_neighbor( cylinder_projections )
%PERIODIC_ Summary of this function goes here
%   Detailed explanation goes here

N = length(cylinder_projections.z_position);

z_pos = cylinder_projections.z_position;
a_pos = cylinder_projections.arc_position;
a_l = cylinder_projections.arclength;

x_periodic = [z_pos, z_pos, z_pos; a_pos, a_pos-a_l, a_pos+a_l]';

if (size(x_periodic,2) ~= 2)
    msg = 'Input matrix has wrong dimensions! Must be Nx2.'; 
    error(msg)
end

y = pdist(x_periodic);
Y = squareform(y);
Y(logical(eye(size(Y)))) = Inf; 
for i=1:N
    Y(N+i,i) = Inf;
    Y(i,N+i) = Inf;
end

[nearest_dist,neighbor_idx] = min(Y);


neighbor_idx = mod(neighbor_idx,N);
neighbor_idx(neighbor_idx==0) = N;
neighbor_idx = neighbor_idx(1:N); 
nearest_dist = nearest_dist(1:N);



end

