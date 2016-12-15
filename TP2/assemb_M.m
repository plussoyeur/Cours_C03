function [M] = assemb_M(mesh)
% M=ASSEMB_M(K,mesh) assemble et retourne la matrice de masse EF P1 de Lagrange 
% sur le maillage mesh ou mesh est une structure contenant les champs 
% nbs,nbt,elm_som,som_coo,som_zon
%
% Ne tient pas compte des Conditions aux Limites
% K est suppose constant par element et transmis sous forme de tableau
% colonne à mesh.nbt lignes
%



% Initialisation
M = sparse(mesh.nbs, mesh.nbs);

% Boucle sur les triangles pour le calcul des AKij
for ie = 1:mesh.nbt
    
    is = mesh.elm_som(ie,:);
    aire = .5*abs(det([mesh.som_coo(is,:),ones(3,1)]));
    
    
    Mlm = aire/3.0*(0.25*ones(3,3)+.25*eye(3,3));    
    M(is,is) = M(is,is) + Mlm; 
    
end;
    