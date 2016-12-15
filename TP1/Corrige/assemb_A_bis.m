function [A] = assemb_A(K, mesh)
% A=ASSEMB_A(K,mesh) assemble et retourne la matrice de rigidite EF P1 de Lagrange 
% c'est a dire la matrice EF associe a l'operateur -grad (K grad U)
% sur le maillage mesh ou mesh est une structure contenant les champs 
% nbs,nbt,elm_som,som_coo,som_zon
%
% Ne tient pas compte des Conditions aux Limites
% K est suppose constant par element et transmis sous forme de tableau
% colonne à mesh.nbt lignes
%
% Allocation memoire d'une matrice creuse
A  = sparse(mesh.nbs,mesh.nbs);

% Assemblage de la matrice
for ie = 1:mesh.nbt
  is = mesh.elm_som(ie,:); % numeros globaux de sommets
  x  = mesh.som_coo(is,:); % coordonnees sommets 1 2 3
  
  mes = 0.5 * (x(2,1)-x(1,1))*(x(3,2)-x(1,2)) ...
      - 0.5 * (x(2,2)-x(1,2))*(x(3,1)-x(1,1));
  
  a  = ( -x([ 3 1 2],2) + x([2 3 1],2) ) /(2*mes);
  b  = (  x([ 3 1 2],1) - x([2 3 1],1) ) /(2*mes);  
  %mes = 0.5 * det([x, [1;1;1]]);
  
  % Matrice elementaire 3x3
  % i,j = 1, 2, 3
  A_elm = mes * K(ie) * (a'*a + b'*b) ; % matrice 3x3 
       % A_{ij}^K = |K| * kappa(K) * (a_i a_j + b_i b_j)
  
  % Contribution matrice globale
  A(is,is) = A(is,is) + A_elm;
end;