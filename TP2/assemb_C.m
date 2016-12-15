function C = assemb_C(beta, mesh)

% Intialisation
C = sparse(mesh.nbs, mesh.nbs);

% Boucle sur les triangles
for ie = 1:mesh.nbt
    
    is = mesh.elm_som(ie,:);
   
    x = mesh.som_coo(is,1);
    y = mesh.som_coo(is,2);
   
    
    aire = 0.5*abs(det([x, y, ones(3,1)]));
    
    a = (1/(2*aire))*(circshift(y,2)-circshift(y,1));
    b = (1/(2*aire))*(circshift(x,1)-circshift(x,2));  
    c = 1-a.*x-b.*y;
    
    xG = sum(x)/3;
    yG = sum(y)/3;
    
    delta_part = (beta(1)*a + beta(2)*b);
    phi_part   = a*xG+b*yG+c;
    
    Clm = aire*(delta_part'*phi_part);
    
    C(is,is) = C(is,is) + Clm;
    
end


end


