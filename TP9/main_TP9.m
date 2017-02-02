clear all;

%-- Lecture du maillage
mesh = lect_mesh('carre11');
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);
tri = mesh.elm_som;

%-- Champ de vitesse tournant
u = @(x,y) [ -sin(pi*y)*cos(.5*pi*x) ; sin(pi*x)*cos(.5*pi*y) ];

%-- Donnee initiale
y0 = zeros(mesh.nbs,1);
yfonc = @(x,y) max( 4*(0.25 - sqrt( (x-.5).^2 + y.^2)  ), 0);
ysol = yfonc(x,y);

%-- Donnees
Gammam = [];
dt = .01;
g = @(x,y) 1;
t = 0;
T = pi;

ysol = convect(mesh,u,ysol,Gammam,g,dt);
trisurf(tri,x,y,ysol);


%-- Iteration
while(t < T)
   t = t + dt; 
   ysol = convect(mesh,u,ysol,Gammam,g,dt);
   
%    clf();
%    subplot(1,2,1), trisurf(tri,x,y,ysol);
%    subplot(1,2,2), tri_contour(tri,x,y,ysol, 0:0.05:1); grid;
%    %trisurf(tri,x,y,ysol)
%    axis([-1 1 , -1 1]);
%    view(-29,46);
%    drawnow();
   
clf();
h_tri = trisurf(tri, x, y, ysol);
light
lighting gouraud
material dull
shading interp
set(h_tri, 'EdgeColor', 'none');
axis([-1 1 -1 1 0 1]);
drawnow();

end