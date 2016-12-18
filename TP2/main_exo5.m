%  Perte de la coercivite

mesh = lect_mesh('square');
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);

% Coefficient de diffusion
kappa = ones(mesh.nbt,1);
beta = ones(2,1);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Test de l'assemblage de C
% mesh2 = lect_mesh('square');
% 
% x = mesh2.som_coo(:,1);
% y = mesh2.som_coo(:,2);
% U = x + y;
% 
% C = assemb_C(beta,mesh2);
% U'*C*U
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
