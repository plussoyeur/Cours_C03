function [ sol ] = init( u0, mesh )
% Calcule la valeur de u0 pour chaque triangle K au centre de gravite du
% triangle et stocke ces resultats dans un vecteur sol de taille le nombre
% de triangle

sol = u0(mesh.elm_gra(:,1),mesh.elm_gra(:,2));

end

