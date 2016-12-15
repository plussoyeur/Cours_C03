%  Resolution du systeme element fini de l'equation de Poisson avec
%  conditions aux limites de Dirichlet

mesh = lect_mesh('DOM1');
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);

% Coefficient de diffusion
kappa = ones(mesh.nbt,1);

% Declaration de f
f = @(x,y) (0.*x+0.0*y); % fonction anonyme

% Declaration de g
% g = zeros(mesh.nbs,1);
% for i = 1:mesh.nbs
%     if mesh.som_zon(i) == 2
%         g(i) = 3;
%     elseif mesh.som_zon(i) == 3
%         g(i) = -1;
%     end
% end
g = @(z,x,y) ((z==1)*(1.0)+(z==2)*(2.0)+(z==3)*(-1.0));

% Assemblage des matrices
A = assemb_A(kappa,mesh);
F = assemb_F(f,mesh);

% Initialisation de u
u = zeros(mesh.nbs,1);

% Assemblage du second membre complet
dir = find(mesh.som_zon ~= 0);
inconnues = setdiff((1:mesh.nbs)',dir);
u(dir) = g(mesh.som_zon(dir),x(dir),y(dir));
F = F - A*u;

% Resolution du systeme
u(inconnues) = A(inconnues, inconnues)\F(inconnues);

% Representation graphique
tri = mesh.elm_som;
figure(1); clf; colormap(cool);
trimesh(tri,x,y,u);

% % Verification de l'assemblage de M avec le maillage square.ambda
% mesh2 = lect_mesh('square');
% M = assemb_M(mesh2);
% V = ones(mesh2.nbs,1);
% V'*M*V



