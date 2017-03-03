function f = voronoi2d( cylinder_projections )
%VORONOI2D Summary of this function goes here
%   Detailed explanation goes here

nPoints = 500; 

zpos = linspace(cylinder_projections.zmin, cylinder_projections.zmax, nPoints);
al = cylinder_projections.arclength;
arcpos = linspace(0, al, nPoints);

area_per_pixel = (zpos(2) - zpos(1))*(arcpos(2) - arcpos(1));
[ZPOS, ARCPOS] = meshgrid(zpos, arcpos);
% size([zpos_periodic; arcpos_periodic]')
size([ZPOS(:) ARCPOS(:)])

arcpos_p = cylinder_projections.arc_position;
zpos_p = cylinder_projections.z_position;

arcpos_periodic = [arcpos_p - al, arcpos_p, arcpos_p + al];
zpos_periodic = [zpos_p, zpos_p, zpos_p];

[IDX, dist] = knnsearch([zpos_periodic; arcpos_periodic]',[ZPOS(:) ARCPOS(:)]);
IDX = reshape(IDX, nPoints, nPoints);
IDX = mod(IDX, length(zpos_p));

dist = reshape(dist, nPoints, nPoints);

f1 = figure;
imagesc(zpos, arcpos, IDX)
set(gca, 'ydir', 'normal')
hold on
plot(zpos_p, arcpos_p, 'wo', 'markerfacecolor', 'w')
xlabel('Longitudinal position (um)')
ylabel('Arc position (um)')
axis equal
title('Voronoi diagram of fiber surface')


histlist = [];
for i=min(min(IDX)):max(max(IDX))
    IDX2 = IDX==i;
    
    [row, col] = find(IDX2);
    if min(col) > 1 && max(col) < nPoints
        s = sum(sum(IDX2));
        if s > 0
            histlist(end+1) = s*area_per_pixel;
        end
    end
end

f2 = figure;
f = boxplot(histlist);
%histogram(histlist, 'Normalization', 'probability')
title('Domain area distribution (um x um)')
%xlabel('Surface area (um x um)')
ylabel('Probability')
histlist;


% IDX2(1,:) = 100; 
% IDX2(end,:) = 100;


f3 = figure;
imagesc(zpos, arcpos, dist)
set(gca, 'ydir', 'normal')
xlabel('Longitudinal position (um)')
ylabel('Arc position (um)')
title('Distance map of fiber surface')
axis equal
colorbar



f4 = figure;
histogram(dist(:), 'Normalization', 'probability')
xlabel('Distance (um)')
ylabel('Probability')
title('Distance distribution for a random point to its nearest nuclei')

f = {};
f.voronoi = f1;
f.areahist = f2; 
f.dist = f3; 
f.disthist = f4; 
f.arealist = histlist;
end

