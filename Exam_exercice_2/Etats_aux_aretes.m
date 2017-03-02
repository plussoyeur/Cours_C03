function [uKL_m, uKL_p] = Etats_aux_aretes( mesh,sol_t,gK,alphaK_KL,alphaL_KL )
%Emeriau Pierre-Emmanuel
%emeriau.pe@gmail.com

alpha = zeros(mesh.nbt);
uKL_m = zeros(mesh.nba);
uKL_p = zeros(mesh.nba);

for it = 1:mesh.nbt
    ia = mesh.elm_fac(it,:);
    aplha(it) = min(alphaK_KL(ia));
end

for ia = 1:mesh.nba
    ie = mesh.fac_elm(ia,:); % liste des numeros de K et L
    
    xK = mesh.elm_gra(ie(1),:);
    xKL = mesh.fac_gra(ia);  
    uK = sol_t(ie(1));

    
    uKL_m = uK + alpha(ie(1))*gK(ie(1),:)*(xKL-xK)';
    
    if(ie(2) ~= 0)
        uL = sol_t(ie(2));
        xL = mesh.elm_gra(ie(2),:);
        uKL_p = uL + alpha(ie(2))*gK(ie(2),:)*(xKL-xL)';
    end
    
end


end

