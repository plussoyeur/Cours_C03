function [ sol_s ] = tri_to_sum( mesh, sol_t )
% En entrée : sol_t (solution aux triangles) -> taille mesh.nbt
% En sortie : sol_s (solution aux sommets) -> taille mesh.nbs

    sol_s = zeros(mesh.nbs,1);
    num = zeros(mesh.nbs,1);
    denom = zeros(mesh.nbs,1);

    for ie = 1:mesh.nbt
       
        is = mesh.elm_som(ie,:);
        
        num(is) = num(is) + mesh.elm_mes(ie)*sol_t(ie)*ones(3,1);
        denom(is) = denom(is) + mesh.elm_mes(ie)*ones(3,1);
        
        
    end
    
    sol_s = num./denom;
end

