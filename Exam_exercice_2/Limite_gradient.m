function [ alphaK_KL, alphaL_KL ] = Limite_gradient( mesh,sol_t,gK )
% Emeriau Pierre-Emmanuel
% emeriau.pe@gmail.com

alphaK_KL = ones(mesh.nba);
alphaL_KL = ones(mesh.nba);

for ia = 1:mesh.nba
    ie = mesh.fac_elm(ia,:); % liste des numeros de K et L
    
    xK = mesh.elm_gra(ie(1),:);
    xKL = mesh.fac_gra(ia);  
    uK = sol_t(ie(1));
    
    if(ie(2) ~= 0)
        xL = mesh.elm_gra(ie(2),:);
        uL = sol_t(ie(2));
    else
        uL = uK;
    end
    
    for j = 1:8
       t = uK + alphaK_KL(ia)*gK(ie(1),:)*(xKL-xK)';
       if(t <= min(uK,uL) || t >= max(uK,uL))
           alphaK_KL(ia) = .5*alphaK_KL(ia);
       else
           break
       end
    end
    
    if(ie(2) ~= 0)
        for j = 1:8
            t = uL + alphaL_KL(ia)*gK(ie(2),:)*(xKL-xL)';
            if(t <= min(uK,uL) || t >= max(uK,uL))
                alphaL_KL(ia) = .5*alphaL_KL(ia);
            else
                break
            end
        end
    end
    
    
    
       
end

end

