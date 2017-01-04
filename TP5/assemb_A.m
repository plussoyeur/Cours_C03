function [A] = assemb_A(K,mesh)
% A=ASSEMB_A(K,mesh) assemble et retourne la matrice de rigidite EF P1 de Lagrange 
% c'est a dire la matrice EF associe a l'operateur -grad (K grad U)
% sur le maillage mesh ou mesh est une structure contenant les champs 
% nbs,nbt,elm_som,som_coo,som_zon
%
% Ne tient pas compte des Conditions aux Limites
% K est suppose constant par element et transmis sous forme de tableau
% colonne à mesh.nbt lignes
%

% Intialisation
A = sparse(mesh.nbs, mesh.nbs);


% Boucle sur les triangles
for ie = 1:mesh.nbt
    
    is = mesh.elm_som(ie,:);
    
    a = circshift(mesh.som_coo(is,2),2)-circshift(mesh.som_coo(is,2),1);
    b = circshift(mesh.som_coo(is,1),1)-circshift(mesh.som_coo(is,1),2);
    
    aire = 0.5*abs(det([mesh.som_coo(is,:), ones(3,1)]));
    
    Alm = (K(ie)/(4*aire))*(a*a' + b*b');
    
    A(is,is) = A(is,is) + Alm;
    
end