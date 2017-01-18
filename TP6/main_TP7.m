clear all;

% Declaration du maillage utilise
mesh = lect_mesh('disq0');
mesh = face_number(mesh); % Structure de donnees VF
%

% Donnes initiales
T = pi;
deltat = 0.015;
u0 = @(x,y)(exp(-50*((x-0.4).^2+(y-0.0).^2)) ); % donnee initiale exp
u0 = @(x,y)( 0+((x-0.4).^2+(y-0.0).^2)<=0.04);
%
snew_t = init(u0,mesh)
%

% Boucle sur le temps
t = 0;
while(t < T)
    t = t + deltat;
    
    snew_t = conv_sca(mesh,snew_t,deltat);
    
    % Verification et visualisation
    snew_s = tri_to_sum(mesh, snew_t);
    tri = mesh.elm_som;
    x = mesh.som_coo(:,1);
    y = mesh.som_coo(:,2);
    clf();
    
    trisurf(tri, x, y, snew_s);
    %subplot(1,2,1), trisurf(tri, x, y, snew_s);
    %subplot(1,2,2), tri_contour(tri, x, y, sol_s, 0:0.05:1); grid;
    drawnow();  
end;