function ynew = convect(mesh, u, y, Gammam, g, dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fonction qui met en oeuvre une methode EF des caracteristiques
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-- Initialisation 
ynew = zeros(mesh.nbs,1);
flag = zeros(mesh.nbs,1);

%-- On attribue les valeurs aux sommets Gammam
ynew(Gammam) = g(Gammam);
flag(Gammam) = 1;

%-- Boucle sur les triangles
for ie = 1:mesh.nbt
    is = mesh.elm_som(ie,:);
    %-- Boucle sur les sommets du triangle
    xi = mesh.som_coo(is,1);
    yi = mesh.som_coo(is,2);    
    S = [mesh.som_coo(is(1),:); mesh.som_coo(is(2),:); mesh.som_coo(is(3),:)];
    A = [[S(1,:)';1],[S(2,:)';1],[S(3,:)';1]];
    yold = [y(is(1)); y(is(2)); y(is(3))];
    %A = [xi(is(1)) xi(is(2)) xi(is(3)); ...
    %     yi(is(1)) yi(is(2)) yi(is(3)); ...
    %     1          1         1];
    for i = 1:3
        if (flag(is(i)) == 1)
            continue;
        elseif (flag(is(i)) == 0)
            % Traitement du sommet
            P_i = S(i,:)' - dt*u( xi(i) , yi(i) );
            % Second membre
            B = [P_i; 1];
            lambda = A \ B;
            if (min(lambda) >= 0)            
                ynew(is(i)) = lambda'*yold;
                flaf(is(i)) = 1;
            end
        else
            disp('ERROR : flag problem. Should either 1 or 0');
        end
    end
end

end

