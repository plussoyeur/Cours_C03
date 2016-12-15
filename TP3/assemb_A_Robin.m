function [A] = assemb_A_Robin(K,alpha,mesh)
% A=ASSEMB_A(K,mesh) assemble et retourne la matrice de rigidite EF P1 de Lagrange 
% c'est a dire la matrice EF associe a l'operateur -grad (K grad U)
% sur le maillage mesh ou mesh est une structure contenant les champs 
% nbs,nbt,elm_som,som_coo,som_zon
%
% Ne tient pas compte des Conditions aux Limites
% K est suppose constant par element et transmis sous forme de tableau
% colonne à mesh.nbt lignes
%
% Modifications pour la prise en compte de CONDITIONS DE ROBIN
%


% Initialisation
A = sparse(mesh.nbs, mesh.nbs);

% Boucle sur les triangles pour le calcul des AKij
for ie = 1:mesh.nbt
    
    is = mesh.elm_som(ie,:);
    
    a = circshift(mesh.som_coo(is,2),2)-circshift(mesh.som_coo(is,2),1);
    b = circshift(mesh.som_coo(is,1),1)-circshift(mesh.som_coo(is,1),2);
    
    aire = .5*abs(det([mesh.som_coo(is,:),ones(3,1)]));
    
    Alm = K(ie)/(4*aire)*(a*a' + b*b');  
    
    id = is(find(mesh.som_zon(is) ~= 0));
    if(max(size(id)) == 2)
        x = mesh.som_coo(id(1),1) - mesh.som_coo(id(2),1);
        y = mesh.som_coo(id(1),2) - mesh.som_coo(id(2),2);
        d = sqrt(x^2+y^2);
        A(id,id) = A(id,id) + alpha*d/3*(.5*ones(2,2)+.5*eye(2,2));
    end
    
    A(is,is) = A(is,is) + Alm; 
    
    
end;
    