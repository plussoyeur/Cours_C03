function out=trace_bord(mesh,col)
% TRACE_BORD(MESH) trace le bord du domaine en bleu
% MESH est une structure dont les champs sont
% nbt,nbs,nbab,elm_som,elm_coo,som_zon,abd_som
%
% TRACE_BORD(MESH,COL) 
% ou col est une chaine compose de l'un des 
% caracteres y,m,c,r,g,b,w,k trace le bord avec la couleur col
%
% Copyright (c) 2005 by Frederic Pascal, ENS de Cachan 

if nargin==2 & ischar(col),
  color = col;
else
  color = 'b';
end

% Recuperation des sommets des aretes
aretes = mesh.abd_som(:,[1 2])';
nba = size(aretes,2);

% Traces des sommets
fig=plot(reshape(mesh.som_coo(aretes,1),2,nba),...
         reshape(mesh.som_coo(aretes,2),2,nba),color);
     
% Adaptation de la fenetre graphique
axis equal
axis tight

if nargout == 1, out = fig; end



