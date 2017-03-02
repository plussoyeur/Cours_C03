% Emeriau Pierre-Emmanuel
% emeriau.pe@gmail.com
clear all;

% Lecture du maillage
mesh = lect_mesh('Disque');
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);
tri = mesh.elm_som;

% choix de kappa
kappa = ones(mesh.nbt,1);

% donnees
dt = 0.01;
T = 3;

%%%%%%%%%%% Calcul de u0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
u0 = zeros(mesh.nbs,1);

% Assemblage du second membre complet
dir = find(mesh.som_zon == 1);
inconnues = setdiff((1:mesh.nbs)',dir);
u0(dir) = g(x(dir),y(dir));
F = F - A*u0;

% Resolution du systeme
u0(inconnues) = A(inconnues, inconnues)\F(inconnues);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Declaration du second membre fchaleur
fchaleur = @(x,y) (0);

% Dirichlet 
udir = zeros(mesh.nbs,1);
udir(dir) = g(x(dir),y(dir));

% Assemblage de la matrice de masse
M = assemb_M(mesh);

% Matrice devant u_(n+1)
A_tilde = dt*A + M;

it = 0;
u = u0;
t=0;
trimesh(tri, x, y, u);

while(t<T)
    it = it+1;
    t = t+dt;
    % Assemblage de f
    F = assemb_F(fchaleur,mesh);
    F = F + M*u;
    
    % Pseudo elimination
    F = F-A_tilde*udir;
    
    % Resolution du système lineaire
    u(inconnues) = A_tilde(inconnues, inconnues)\F(inconnues);
     
    if(mod(it,10) == 0)  
        figure(1); clf; colormap(cool);
        trimesh(tri, x, y, u);
        view(-16,8); 
        drawnow();
        pause(1);
    end
end 





