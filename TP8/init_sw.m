function [Vnew] = init_sw(mesh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init_sw initialise l'inconnue V
% V(:,1) : hauteur d'eau
% V(:,2) : hauteur d'eau x vitesse en x
% V(:,3) : hateur d'eau x vitesse en y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Definition des vecteurs x et y du maillage
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);

% Initialisation de V
Vnew = zeros(3,mesh.nbt);

% Liste des triangles de part et d'autres du barrage
x_som_1 = x(mesh.elm_som(:,1)); % coordonnes du 1er sommet dans un triangle
x_som_2 = x(mesh.elm_som(:,2)); % ...
x_som_3 = x(mesh.elm_som(:,3)); % ...

up_bar = find(max(max(x_som_1,x_som_2),x_som_3)<=100);
down_bar = setdiff(1:mesh.nbt,up_bar);

% Hauteur d'eau initiale de part de d'autre du barrage
Vnew(1,up_bar) = 10;
Vnew(1,down_bar) = 5;


end

