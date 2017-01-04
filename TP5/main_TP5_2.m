clear all;

% Lecture du mesh
mesh = raf_mesh(lect_mesh('L2'));

% choix de kappa
kappa = 0.01*ones(mesh.nbt,1);

% choix de tau
tau = 1;

% choix de T
T = 5;

% choix de delta_t
delta_t = 0.1;

% définition de la fonction f
f = @(x,y,t) sin(2*pi*(x+y-t/tau));

% Assemblage de la matrice de rigidité
A = assemb_A(kappa, mesh);

% Assemblage de la matrice de masse
M = assemb_M(mesh);

    
% Matrice devant u_(n+1)
A_tilde = A + (1/delta_t)*M;


% Recuperation des donnees au bord
dir = find(mesh.som_zon == 1);
inconnues = setdiff(1:mesh.nbs, dir);
u(dir) = 0;

% Données itérations max
maxIter = floor(T/delta_t);

% Initialisation de l'inconnue
u = zeros(mesh.nbs,1);

% Création de la figure
h_fig = figure(1);

for n = 0:maxIter

    t = delta_t*n;
    
    % Assemblage de f
    F = assemb_F_time(f,mesh,t);
    F = F + 1/delta_t*M*u;
    
    % Pseudo elimination
    F = F-A_tilde*u;
    
    % Resolution du système lineaire
    u(inconnues) = A_tilde(inconnues, inconnues)\F(inconnues);
     
    tri = mesh.elm_som;
    x = mesh.som_coo(:,1);
    y = mesh.som_coo(:,2);
    
    pause(0.1);    
    trimesh(tri, x, y, u);
    view(4.5,70);    
end 