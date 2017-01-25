function [ A ] = matrice_A( nu, V )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fonction qui retourne la matrice A d'interet pour calculer le signe de M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Definition des variables
u = zeros(2,1);
u = V(2:3)/V(1);
h = V(1);

g = 9.81;
c = sqrt(g*h);

% Calcul de A
A = [0, nu(1), nu(2) ; ...
    (c^2 - u(1)^2)*nu(1) - u(1)*u(2)*nu(2), 2*u(1)*nu(1)+u(2)*nu(2), u(1)*nu(2); ...
    (c^2 - u(2)^2)*nu(2) - u(2)*u(1)*nu(1), u(2)*nu(1), 2*u(2)*nu(2)+u(1)*nu(1)
    ];


end

