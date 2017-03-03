function plotline_hObject_thing(hObject,event,ellipse, X,Y,Z,data1, data2,plot_pointer, index_list)
n = get(hObject,'Value');
n = max(round(n*length(index_list)),1);
n = max(max(index_list)-n+1, 1);
clim1 = [min(data1(:)), max(data1(:))];
clim2 = [min(data2(:)), max(data2(:))];
zlist = linspace(min(Z(:)), max(Z(:)),10);
philist = linspace(0,2*pi,100); 
xlist = ellipse.a*cos(philist); 
ylist = ellipse.b*sin(philist);
xvec = [xlist; ylist]; 

% S = [cos(ellipse.phi), -sin(ellipse.phi); sin(ellipse.phi), cos(ellipse.phi)];
% xvec = S*xvec; 


% xmin = min(X(:)); ymin = min(Y(:)); zmin = min(Z(:));
% xmax = max(X(:)); ymax = max(Y(:)); zmax = max(Z(:));

X = X(:,:,1:n);
Y = Y(:,:,1:n);
Z = Z(:,:,1:n);
data1 = data1(:,:,1:n);
data2 = data2(:,:,1:n);

xmin = min(X(:)); ymin = min(Y(:)); zmin = min(Z(:));
xmax = max(X(:)); ymax = max(Y(:)); zmax = max(Z(:));

nPoints3d = length(X);

[Xe,Ye,Ze] = cylinder(linspace(1,1,nPoints3d),300);
Ze = Ze*(zmax-zmin) + zmin;
Xe = Xe*ellipse.a;
Ye = Ye*ellipse.b; 

% S = [cos(ellipse.phi), -sin(ellipse.phi); sin(ellipse.phi), cos(ellipse.phi)];
% 
% for i=1:length(Xe(:,1))
%     for j=1:length(Xe(1,:))
%         xvect = [Xe(i,j); Ye(i,j)]; 
%         xvect = S*xvect;
%         Xe(i,j) = xvect(1);
%         Ye(i,j) = xvect(2);
%     end
% end

[Xt,Yt,Zt] = sphere(nPoints3d);
Zt(:,:) = zmax;
Xt = Xt*ellipse.a;
Yt = Yt*ellipse.b; 
 
[Xb,Yb,Zb] = sphere(nPoints3d);
Zb(:,:) = zmin;
Xb = Xb*ellipse.a;
Yb = Yb*ellipse.b; 

figure(plot_pointer)



subplot(1,2,1)
hold off
slice(X,Y,Z,data1, Xe, Ye, Ze);

hold on

slice(X,Y, Z, data1,Xt, Yt, Zt )

slice(X,Y, Z, data1,Xb, Yb, Zb )

a = get(gca, 'children');
for i=1:length(a)
   set(a(i), 'edgecolor', 'none')
end
for i=1:length(zlist)
   plot3(xvec(1,:), xvec(2,:), linspace(zlist(i), zlist(i), 100), 'k')
end
caxis(clim1)

camlight headlight
lighting phong
title('3D Voronoi diagram')
xlabel('x [um]')
ylabel('y [um]')
zlabel('z [um]')
subplot(1,2,2)
hold off
slice(X,Y,Z,data2, Xe, Ye, Ze);

hold on

slice(X,Y, Z, data2,Xt, Yt, Zt )

slice(X,Y, Z, data2,Xb, Yb, Zb )

a = get(gca, 'children');
for i=1:length(a)
   set(a(i), 'edgecolor', 'none')
end
for i=1:length(zlist)
   plot3(xvec(1,:), xvec(2,:), linspace(zlist(i), zlist(i), 100), 'k')
end
caxis(clim2)
camlight headlight
lighting phong
title('3D distance map')
xlabel('x [um]')
ylabel('y [um]')
zlabel('z [um]')

end

