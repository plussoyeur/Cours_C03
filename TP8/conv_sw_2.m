function [Vnew_t] = conv_sw_2(mesh, Vold_t, dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fonction qui returne la solution a l'itere (n+1) du schema a partir de
%     l'itere au temps n
% En entree : ancien ietere Vold
%             maillage
%             pas de temps dt
% En sortie : nouvel itere Vnew
%
% Nouveaute : on utilise le schema de Lagrange-Flux au lieu de l'ancien
%             flux            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

% -- Grandeur du probleme 
g = 9.81;

% -- Normal aux aretes
N1_KL = mesh.fac_nor(:,1);
N2_KL = mesh.fac_nor(:,2);

% -- Initialisation
Vnew_t = zeros(3,mesh.nbt);
dsol_t = zeros(3,mesh.nbt);

% -- Boucle sur les aretes
for ia=1:mesh.nba
    
    % -- Element K et L sur l'autre ia
    ie = mesh.fac_elm(ia,:);
    
    nu = [N1_KL(ia); N2_KL(ia)];
    
    % Mesure de l'arete
    mesA = mesh.fac_mes(ia);
    
    % -- Flux pour aretes interieures
    if(mesh.fac_zon(ia) == 0)
        % Hauteur pour les triangles K et L
        hK = Vold_t(1,ie(1));
        hL = Vold_t(1,ie(2));
        % Celerite par triangle K ou L
        cK = sqrt(g*hK);
        cL = sqrt(g*hL);
        % Pression par triangle K ou L
        pK = .5*g*hK^2;
        pL = .5*g*hL^2;
        % Champs aux triangles consideres
        U_K = Vold_t(:,ie(1));
        U_L = Vold_t(:,ie(2));
        
        % Calcul des elements propres au schema Lagrange-Flux
        p_star = (pK*hL+pL*hK)/(hL+hK) - ...
            hK*hL*max(cK,cL)/(hL+hK) * ...
            ( U_L([2:3])/hL - U_K([2:3])/hK )'*nu;
        
        unu_star = .5*( U_K([2:3])/hK + U_L([2:3])/hL )'*nu - ...
            1.0/((hL+hK)*max(cK,cL)) * (pL-pK);
        
        phi =  U_K*max(0,unu_star) + U_L*min(0,unu_star) + ...
            [0; p_star*nu(1); p_star*nu(2) ];

               
        dsol_t(:,ie(1)) = dsol_t(:,ie(1)) + mesA*phi;
        dsol_t(:,ie(2)) = dsol_t(:,ie(2)) - mesA*phi;
    end
    
    % -- Flux pour aretes de bord de type 1
    if(mesh.fac_zon(ia) == 1)
        phi =  0.5*[0; ...
                    g*Vold_t(1,ie(1))^2*nu(1); ...
                    g*Vold_t(1,ie(1))^2*nu(2)];
        dsol_t(:,ie(1)) = dsol_t(:,ie(1)) + mesA*phi;        
    end
    
    % -- Flux pour arete de bord artificiel de type 2
    if(mesh.fac_zon(ia) == 2)
        ps = Vold_t(2:3,ie(1))'*nu;
        h = Vold_t(1,ie(1));
        phi =  [ps ; ...
                ps*Vold_t(2,ie(1))/h + .5*g*h^2*nu(1) ; ...
                ps*Vold_t(3,ie(1))/h + .5*g*h^2*nu(2) ;
               ];
        dsol_t(:,ie(1)) = dsol_t(:,ie(1)) + mesA*phi;        
    end
    
    
end

% -- Mise a jour du resultat

mes = [mesh.elm_mes'; mesh.elm_mes'; mesh.elm_mes'];

Vnew_t = Vold_t - dt*dsol_t./mes;

end
