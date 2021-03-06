function [C] = assemb_C(beta,mesh)
% Initialisation
C = sparse(mesh.nbs, mesh.nbs);

% Boucle sur les triangles pour le calcul des AKij
for ie = 1:mesh.nbt
    
    is = mesh.elm_som(ie,:);
    x = mesh.som_coo(is,1);
    y = mesh.som_coo(is,2);
    
    aire = .5*abs(det([x,y,ones(3,1)]));
    
    a = (1/(2*aire))*(circshift(y,2)-circshift(y,1));
    b = (1/(2*aire))*(circshift(x,1)-circshift(x,2));
    
    xG = sum(x,1)/3;
    yG = sum(y,1)/3;
    
    Clm = aire/3*(beta(1)*a+beta(2)*b)*ones(1,3);
    
    C(is,is) = C(is,is) + Clm; 
    
   
end;
    