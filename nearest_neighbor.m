function [ neighbor_idx, nearest_dist, dist, f ] = nearest_neighbor( cylinder_projections )
%NEAREST_NEIGHBOR(x) finds the distances and indices of nearest neighbors
%in a set of points, as well as a matrix containing all distances
%   Author: Andreas Våvang Solbrå
%   Date: 01-12-2015 (dd-mm-yyyy)
%
%   Description: 
%   Input matrix must be Nx3
%
%   Example: 
%   x = rand(4,3); 
%   [idx,nearest_dist, dist] = nearest_neighbor(x);
%   
x = [cylinder_projections.z_position; cylinder_projections.arc_position]';

if (size(x,2) ~= 2)
    msg = 'Input matrix has wrong dimensions! Must be Nx2.'; 
    error(msg)
end

y = pdist(x);
Y = squareform(y);
Y(logical(eye(size(Y)))) = Inf; 
[nearest_dist,neighbor_idx] = min(Y);
Y(logical(eye(size(Y)))) = 0;
dist=Y;

f = 0;
% 
% figure
% f = boxplot(nearest_dist);
% 
% title('Nearest neighbor distance distribution')
% xlabel('Distance (um)')
% ylabel('Probability')
end

