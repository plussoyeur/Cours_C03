function [ gK ] = Calcule_gradient(mesh,sol_t)
% Emeriau Pierre-Emmanuel
% emeriau.pe@gmail.com

% Initialisation de gk
gK = zeros(mesh.nbt,2);

% Boucle sur les aretes
for ia = 1:mesh.nba
    ie = mesh.fac_elm(ia,:); % liste des numeros de K et L
    nu = mesh.fac_nor(ia,:); % normal a l'arete
          
    % Mesure du triangkle
    mesK1 = mesh.elm_mes(ie(1));
    % Mesure de l'arete
    mesA = mesh.fac_mes(ia);
    % Valeur au milieu de l'arete
    if(ie(2)~=0)
        mesK2 = mesh.elm_mes(ie(2));
        uKL = 0.5*(sol_t(ie(1)) + sol_t(ie(2)) );
        gK(ie(1),:) = gK(ie(1),:) + mesA/mesK1*uKL*nu;
        gK(ie(2),:) = gK(ie(2),:) - mesA/mesK2*uKL*nu;
    else
        uKL = sol_t(ie(1));
        gK(ie(1),:) = gK(ie(1),:) + mesA/mesK1*uKL*nu;
    end


end

