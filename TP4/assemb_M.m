function M = assemb_M(mesh)

% Intialisation
M = sparse(mesh.nbs, mesh.nbs);

% Boucle sur les triangles
for ie = 1:mesh.nbt
    
    is = mesh.elm_som(ie,:);
    
    aire = 0.5*abs(det([mesh.som_coo(is,:), ones(3,1)]));
    
    Mlm = 0.25*ones(3,3)+ 0.25*eye(3,3); 
    Mlm = (aire/3.0)*Mlm;
    
    M(is,is) = M(is,is) + Mlm;
    
end


end

