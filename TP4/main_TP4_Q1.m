clear all;
nb = 20;

% Lecture du mesh
mesh = raf_mesh(raf_mesh(raf_mesh(lect_mesh('L0'))));

% choix de kappa
kappa = ones(mesh.nbt,1);

% Assemblage de la matrice de rigidit�
A = assemb_A(kappa, mesh);
spy(A);

% Assemblage du second membre
M = assemb_M(mesh);

% Initialisation de l'inconnue
u = zeros(mesh.nbs,1);

% Recuperation des donnees au bord
dir = find(mesh.som_zon == 1);
inconnues = setdiff(1:mesh.nbs, dir);

% Elimination des noeuds au bord
A_tronc = A(inconnues,inconnues);
M_tronc = M(inconnues,inconnues);


% Recuperation des vecteurs propres et valeurs propres
[EigVect,EigVal] = eigs(A_tronc,M_tronc,nb,9.0)

phi = zeros(mesh.nbs,6);
[maxVal,maxInd] = max(EigVect);

for k= 1:nb
    phi(inconnues,k) = sign(EigVect(maxInd(k),k))*EigVect(:,k)/maxVal(k);
end

tri = mesh.elm_som;
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);

for j=1:nb
   subplot(4,5,21-j)
   trimesh(tri,x,y,phi(:,j))
end

