function [F] = assemb_F_Robin(fun,alpha,ua,g,mesh)
% F=ASSEMB_F(fun,mesh) assemble le second membre du pb du Laplacien
% avec conditions aux limites de Neumann Homogenes
% c'est-a-dire F_i = \int_{\Omega} fun(x,y)*\phi_i dxdy 
% pour la methode des EF P1 de Lagrange
%
% La fonction fun(x,y) est une fonction de 2 arguments 
%
% Copyright (c) 2005 by Frederic Pascal, 2015 Florian De Vuyst ENS de Cachan

% Allocation
F = zeros(mesh.nbs,1);

% Assemblage du second membre avec une formule a 1 point
for ie = 1:mesh.nbt
  is = mesh.elm_som(ie,:);            % vecteur 1x3
  x  = mesh.som_coo(is,:);            % vecteur 3x2
  a  = -x([ 3 1 2],2)+x([2 3 1],2);   % vecteur 3x1
  b  =  x([ 3 1 2],1)-x([2 3 1],1);
  mesK = 0.5*(b(2)*a(1)-a(2)*b(1));
  
  grav = sum(x,1)/3; % Centre de gravite
  
  F_elm = mesK * fun(grav(1),grav(2)) / 3.0 * ones(3,1);
  
  
  id = is(find(mesh.som_zon(is) ~= 0));
  if(max(size(id)) == 2)
      z = mesh.som_zon(id(1));
      x1 = mesh.som_coo(id(1),1); 
      x2 = mesh.som_coo(id(2),1);
      xm = (x1+x2)/2;
      y1 = mesh.som_coo(id(1),2);
      y2 = mesh.som_coo(id(2),2);
      ym = (y1+y2)/2;
      d = sqrt((x1-x2)^2+(y1-y2)^2);
      F(id) = F(id) + .5*d*(alpha*ua(z,xm,ym)+g(z,xm,ym))*ones(2,1);
  end
    
  
  F(is) = F(is) + F_elm;
end  
