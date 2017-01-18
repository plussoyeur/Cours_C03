function [ snew_t ] = conv_sca( mesh, sold_t, deltat )
% Schema VF upwind explicite d'ordre 1
% En entree : mesh -> maillage VF
%             sold_t -> solution au pas de temps t aux triangles
%             deltat -> pas de temps
% En sortie : snew_t -> solution au pas de temps n+1 aux triangles

snew_t = zeros(mesh.nbt,1);
dsol_t = zeros(mesh.nbt,1);

for ia = 1:mesh.nba % boucle sur les aretes
    ie = mesh.fac_elm(ia,:); % liste des numeros de K et L
    
    if(ie(2) == 0) % on sort de l'iteration si l'arete est sur le bord
       continue 
    end
    
    % vitesse au milieu de l'arete
    c = [-mesh.fac_gra(ia,2),mesh.fac_gra(ia,1)];
    % normal au milieu de l'arete
    nu = mesh.fac_nor(ia,:)*mesh.fac_mes(ia);
    % produit scalaire entre celerite et normal
    ps = c*nu';
    
    if(ps >=0) 
       dsol_t(ie(1)) = dsol_t(ie(1)) + sold_t(ie(1))*ps;
       dsol_t(ie(2)) = dsol_t(ie(2)) - sold_t(ie(1))*ps;
    end

    if(ps <0) 
       dsol_t(ie(1)) = dsol_t(ie(1)) + sold_t(ie(2))*ps;
       dsol_t(ie(2)) = dsol_t(ie(2)) - sold_t(ie(2))*ps;
    end

%     phi = max(0, ps)*sold_t(ie(1)) ...
%         + min(0, ps)*sold_t(ie(2));
%     
%     dsol_t(ie(1)) = dsol_t(ie(1)) + phi;
%     dsol_t(ie(2)) = dsol_t(ie(2)) - phi;
    
end

snew_t = sold_t - deltat*dsol_t./mesh.elm_mes;
