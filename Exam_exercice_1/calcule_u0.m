% Emeriau Pierre-Emmanuel
% emeriau.pe@gmail.com
clear all;

% Lecture du maillage
mesh = lect_mesh('Disque');
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);
tri = mesh.elm_som;

% Declaration du second membre f
f = @(x,y) (20);

% Kappa 
kappa = ones(mesh.nbt,1);

% Dirichlet 
g = @(x,y) (sin(pi.*x)+sin(pi.*y));

% Assemblage des matrices
A = assemb_A(kappa,mesh);
F = assemb_F(f,mesh);

% Initialisation de u
u = zeros(mesh.nbs,1);

% Assemblage du second membre complet
dir = find(mesh.som_zon == 1);
inconnues = setdiff((1:mesh.nbs)',dir);
u(dir) = g(x(dir),y(dir));
F = F - A*u;

% Resolution du systeme
u(inconnues) = A(inconnues, inconnues)\F(inconnues);

% Representation graphique
figure(1); clf; colormap(cool);
trisurf(tri,x,y,u);

