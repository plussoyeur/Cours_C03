clear all;

% Lecture du maillage
mesh = lect_mesh('dam2');
mesh = face_number(mesh); % structure volume fini

% Coordonnees des sommets du maillage
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);
tri = mesh.elm_som;

% Initialisation du vecteur solution
V_t = init_sw(mesh); % solution par triangle

% Choix de dt et de T
T = 6.8;
dt = 0.05;

t = 0;
i = 0;
while(t < T)
    i = i +1;
    t = t + dt;
    V_t = conv_sw(mesh, V_t, dt);
    V_s = tri_to_sum(mesh,V_t(1,:));
    
    % Affichage
    fprintf('Iteration n %i  |  Temps %f  |  \n',i,t);
    fprintf('Hauteur max %f  |  Max vitesse x %f  |  Max vitesse y %f    \n'...
        ,max(V_t(1,:)),max(V_t(2,:)), max(V_t(3,:)));
     fprintf('Hauteur min %f   |  Min vitesse x %f  |  Min vitesse y %f    \n'...
         ,min(V_t(1,:)),min(V_t(2,:)), min(V_t(3,:)))
    
    clf();
    h_tri = trisurf(tri,x,y, V_s);
    light
    lighting gouraud 
    material dull
    %shading interpret
    set(h_tri, 'EdgeColor', 'none');
    view(40,26);
    axis([0,200, 0, 200, 5, 10]);
    drawnow();
end