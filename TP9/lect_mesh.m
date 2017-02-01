function [mesh]=lect_mesh(text)
% MESH=LECT_MESH(TEXT) est une structure associee à un maillage 
% lu dans le fichier dont le nom est la chaine de caracteres [TEXT,'.amdba']  
% Fichier issu du programme EMC2 puis suppression des caracteres -- nbs,nbt
%
% TEST est une chaine de caracteres correspondant au nom du fichier sans
% l'extension
%
% La structure MESH est defini par les champs 
% nbt,nbs,nbab,elm_som,elm_coo,som_zon,abd_som
% 
% Le fichier est de la forme (maillage compose de triangles)
% L_1             nbt  , nbs  
% L_i             i    , x_i  , y_i  , z_i  , zon_i           %coord et zone de sommet i
% L_i+nbs+1       i     ,s_1^i, s_2^i, s_3^i, zon_i           %sommets et zone de triangle i
%
% Copyright (c) 2005 by Frederic Pascal, ENS de Cachan

fid=fopen([text,'.amdba']);
if (fid == -1) 
    error([' LE FICHIER ',text,'.amdba  EST INCONNU '])
end

% Lecture du fichier
aux=dlmread([text,'.amdba']);

% Nb de sommets et de triangles
mesh.nbs  = aux(1,1);
mesh.nbt  = aux(1,2);

% Topologie du maillage
mesh.som_coo = aux(1+1:1+mesh.nbs,2:3);
mesh.som_zon = aux(1+1:1+mesh.nbs,4);
mesh.elm_som = aux(mesh.nbs+1+1:mesh.nbs+1+mesh.nbt,2:3+1);

% calcul et numerotation des aretes de bord
next = [2 3 1];

mesh.nbab = 0;
for ie = 1:mesh.nbt
    for i = 1:3
        is1=mesh.elm_som(ie,i);
        is2=mesh.elm_som(ie,next(i));
        iz1=mesh.som_zon(is1);
        iz2=mesh.som_zon(is2);
        if (iz1~=0) & (iz2~=0) 
            mesh.nbab = mesh.nbab + 1;
            mesh.abd_som(mesh.nbab,1:2) = [is1 is2];
        end
    end
end
