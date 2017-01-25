function [Vnew_t] = conv_sw(mesh, Vold_t, dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fonction qui returne la solution a l'itere (n+1) du schema a partir de
%     l'itere au temps n
% En entree : ancien ietere Vold
%             maillage
%             pas de temps dt
% En sortie : nouvel itere Vnew
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

% -- Grandeur du probleme 
g = 9.81;

% -- Normal aux aretes
N1_KL = mesh.fac_mes(:).*mesh.fac_nor(:,1);
N2_KL = mesh.fac_mes(:).*mesh.fac_nor(:,2);

% -- Initialisation
Vnew_t = zeros(3,mesh.nbt);
dsol_t = zeros(3,mesh.nbt);

% -- Boucle sur les aretes
for ia=1:mesh.nba
    
    % -- Element K et L sur l'autre ia
    ie = mesh.fac_elm(ia,:);
    
    nu = [N1_KL(ia); N2_KL(ia)];
    
    % -- Flux pour aretes interieures
    if(mesh.fac_zon(ia) == 0)
        phi = 0.5*( F_dot_nu(nu, Vold_t(:,ie(1))) + ...                 
                    F_dot_nu(nu, Vold_t(:,ie(2))));
          
        M = matrice_A(nu, .5*(  Vold_t(:,ie(1)) + Vold_t(:,ie(2)) ) );       
        [V,D] = eig(M); 
        SM = V*sign(D) / V; % signe de la matrice M
        
        phi = phi + SM*0.5*( F_dot_nu(nu, Vold_t(:,ie(1))) - ...                 
                    F_dot_nu(nu, Vold_t(:,ie(2))));
    
                
        dsol_t(:,ie(1)) = dsol_t(:,ie(1)) + phi;
        dsol_t(:,ie(2)) = dsol_t(:,ie(2)) - phi;
    end
    
    % -- Flux pour aretes de bord de type 1
    if(mesh.fac_zon(ia) == 1)
        phi =  0.5*[0; ...
                    g*Vold_t(1,ie(1))^2*nu(1); ...
                    g*Vold_t(1,ie(1))^2*nu(2)];
        dsol_t(:,ie(1)) = dsol_t(:,ie(1)) + phi;        
    end
    
    % -- Flux pour arete de bord artificiel de type 2
    if(mesh.fac_zon(ia) == 2)
        ps = Vold_t(2:3,ie(1))'*nu;
        h = Vold_t(1,ie(1));
        phi =  [ps ; ...
                ps*Vold_t(2,ie(1))/h + .5*g*h^2*nu(1) ; ...
                ps*Vold_t(3,ie(1))/h + .5*g*h^2*nu(2) ;
               ];
        dsol_t(:,ie(1)) = dsol_t(:,ie(1)) + phi;        
    end
    
    
end

% -- Mise a jour du resultat

mes = [mesh.elm_mes'; mesh.elm_mes'; mesh.elm_mes'];

Vnew_t = Vold_t - dt*dsol_t./mes;

end
