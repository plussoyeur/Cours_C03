clear all;
%
mesh = lect_mesh('disq0');
mesh = raf_mesh(mesh);
mesh = face_number(mesh); % structure de donnees VF
%
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);
tri = mesh.elm_som;
%
f = @(x,y)(x.*y); % Fonction anonyme
% Calcul de f aux centres des triangles du maillage
sol_t = f(mesh.elm_gra(:,1),mesh.elm_gra(:,2));
sol_s = tri_to_sum(mesh,sol_t);
%
% Trace de la fonction interpolee aux sommets
figure(1);
trisurf(tri,x,y,sol_s);
% Trace de la difference avec solution relle
figure(2);
trisurf(tri,x,y,sol_s-f(mesh.som_coo(:,1),mesh.som_coo(:,2)));


%triplot(mesh.elm_som,x,y);
%trisurf(tri,x,y, f(x,y));

