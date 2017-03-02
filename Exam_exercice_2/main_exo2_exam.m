% Emeriau Pierre-Emmanuel
% emeriau.pe@gmail.com
clear all;

% Lecture du maillage
mesh = lect_mesh('disq1');
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);
tri = mesh.elm_som;

% Enrichissement structures de donnees VF
mesh = face_number(mesh);

% On trace le maillage
trimesh(tri,x,y);

% Verification 1 du calcul du gradient
sol_t = ones(mesh.nbt); % On prend une solution constante
g0 = zeros(mesh.nbt,2);
g0 = Calcule_gradient(mesh,sol_t); % Le gradient doit etre nul
% On trouve bien 0 pour le résultat

% Verifification 2 du calcul du gradient
sol_t = mesh.elm_gra(:,1) + mesh.elm_gra(:,2); % On prend x+y
g1 = zeros(mesh.nbt,2);
g1 = Calcule_gradient(mesh,sol_t); % Le gradient doit valoir 1
% On trouve environ 1 (aux approximations pres)

[alphaK_KL,alphaL_KL] = Limite_gradient(mesh,sol_t,g1);

[uKL_m,uKL_p] = Etats_aux_aretes(mesh, sol_t, g1, alphaK_KL,alphaL_KL);


u0 = init(mesh);
usol = tri_to_sum(mesh,u0);
trisurf(tri,x,y,usol);




