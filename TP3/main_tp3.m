%  Resolution du systeme element fini de l'equation de Poisson avec
%  conditions aux limites de Robin
clear all;

mesh = lect_mesh('DOM2');
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);

% Coefficient de diffusion
kappa = ones(mesh.nbt,1);

% Definition des donnes
alpha = 1;

% Declaration de f
f = @(x,y) (-1); % fonction anonyme

% Declaration de g
g = @(z,x,y) (0);

% Declaration de ua
ua = @(z,x,y)((z==1)*0+(z==2)*2+(z==3)*(-1));

% Assemblage des matrices
A = assemb_A_Robin(kappa,alpha,mesh);
F = assemb_F_Robin(f,alpha,ua,g,mesh);

% Verification de l'assemblage de A avec le maillage square.ambda
%mesh2 = lect_mesh('square');
% A = assemb_A_Robin(kappa,1,mesh);
% V = ones(mesh.nbs,1);
% V'*A*V
% tri = mesh.elm_som;
% trimesh(tri,x,y);


% Initialisation de u
u = zeros(mesh.nbs,1);
% 
% % Assemblage du second membre complet
% dir = find(mesh.som_zon ~= 0);
% inconnues = setdiff((1:mesh.nbs)',dir);
% u(dir) = g(mesh.som_zon(dir),x(dir),y(dir));
% F = F - A*u;
% 
% Resolution du systeme
u = A\F;
% 
% Representation graphique
tri = mesh.elm_som;
figure(1); clf;
trimesh(tri,x,y,u);



