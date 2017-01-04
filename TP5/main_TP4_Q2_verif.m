clear all;

% Lecture du mesh
mesh = raf_mesh(lect_mesh('L0'));

% choix de kappa
kappa = ones(mesh.nbt,1);

% Assemblage de la matrice de rigidité
A = assemb_A(kappa, mesh);

% Assemblage du second membre
M = assemb_M(mesh);
F = assemb_F(@(x,y) 1, mesh); 

% Initialisation de l'inconnue
u = zeros(mesh.nbs,1);
us = zeros(mesh.nbs,1);

% Recuperation des donnees au bord
dir = find(mesh.som_zon == 1);
inconnues = setdiff(1:mesh.nbs, dir);
u(dir) = 0;

% Pseudo elimination
F = F-A*u;

% Resolution du système lineaire
u(inconnues) = A(inconnues, inconnues)\F(inconnues); 


% Elimination des noeuds au bord
A_tronc = A(inconnues,inconnues);
M_tronc = M(inconnues,inconnues);

% Recuperation des vecteurs propres et valeurs propres
[EigVect,EigVal] = eigs(A_tronc,M_tronc,20,9.0);
nb_degen = find(diff(diag(EigVal)) == 0);
if (isempty(nb_degen))
    disp(['Il n''y a pas de valeurs propres degenerees.']);
end
    
phi = zeros(mesh.nbs,20);
[maxVal,maxInd] = max(abs(EigVect));

for k= 1:20
    phi(inconnues,k) = sign(EigVect(maxInd(k),k))*EigVect(:,k)/maxVal(k);
end

A_tilde = phi(inconnues,:)'*A_tronc*phi(inconnues,:);
F_tilde = phi(inconnues,:)'*F(inconnues);

alpha = A_tilde\F_tilde;

for k = 1:20
    us(inconnues) = us(inconnues) + alpha(k).*phi(inconnues,k);
end

tri = mesh.elm_som;
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);

figure(2)
trimesh(tri, x, y, u-us);
