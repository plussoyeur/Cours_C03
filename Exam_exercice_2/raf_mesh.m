function [newmesh]=raf_mesh(oldmesh)
% MESH=RAF_MESH(MESH) raffine le maillage MESH compose de triangles
% en coupant chaque triangle en 4 sous triangles
%
% LES CHAMPS DE MESH DOIVENT ETRE LES SUIVANTS
% mesh.nbs   : nb de sommets
% mesh.nbt   : nb de triangles
% mesh.nbab  : nb d'aretes de bord
%
% mesh.som_coo  : coord des sommets
% mesh.som_zon  : zone des sommets
% mesh.elm_som  : sommets des elements
% mesh.abd_som  : sommets des aretes de bord
%
% Remarque : le numero de zone d'un nouveau sommet = minimun des
%            numeros de zone des 2 sommets voisins de l'arete
%
%Copyright (c) 2005 by Frederic Pascal, ENS de CACHAN
%Copyright (c) 2009 by Jan Valdam, http://sites.google.com/site/janvaldman/software 

% Verification
if isfield(oldmesh,'nbs') ~= 1 ...
        |  isfield(oldmesh,'nbt') ~= 1 ...
        |  isfield(oldmesh,'nbab') ~= 1 ...
        |  isfield(oldmesh,'som_coo') ~= 1 ...
        |  isfield(oldmesh,'som_zon') ~= 1 ...
        |  isfield(oldmesh,'elm_som') ~= 1 ...
        |  isfield(oldmesh,'abd_som') ~= 1
    error([' LES CHAMPS DE LA STRUCTURE MESH NE SONT PAS CORRECT :' ...
            ' RENOMMEZ VOS CHAMPS '])
end

nb_som_elm = size(oldmesh.elm_som,2);
if (nb_som_elm ~= 3) 
    error([' NB SOMMETS PAR ELEMENT NE VAUT PAS 3 '])
end

if isfield(oldmesh,'elm_fac') ~= 1 ...
        | isfield(oldmesh,'fac_som') ~= 1 ...
        | isfield(oldmesh,'nba') ~= 1 
       
    % Faut commencer par numeroter les faces
    % Appel aux programmes de Jan Valdam 
    [oldmesh.elm_fac,oldmesh.fac_som,oldmesh.fac_elm]=getEdges(oldmesh.elm_som);

    % Nba
    oldmesh.nba = size(oldmesh.fac_som,1);
end

% On passe au raffinement
newmesh.nbt = 4*oldmesh.nbt;
newmesh.elm_som = zeros(newmesh.nbt,nb_som_elm);
newmesh.nbs = oldmesh.nbs + oldmesh.nba;
newmesh.som_coo = zeros(newmesh.nbs,2);
if size(oldmesh.som_zon,1) == 1
    newmesh.som_zon = zeros(1,newmesh.nbs);
else
    newmesh.som_zon = zeros(newmesh.nbs,1);
end

% Sommets par elements
newmesh.elm_som(4*(0:oldmesh.nbt-1)+1,:) = oldmesh.elm_fac(1:oldmesh.nbt,:)+oldmesh.nbs;
newmesh.elm_som(4*(0:oldmesh.nbt-1)+2,:) = [oldmesh.elm_fac(1:oldmesh.nbt,1)+oldmesh.nbs ...
    oldmesh.elm_som(1:oldmesh.nbt,3) ...
    oldmesh.elm_fac(1:oldmesh.nbt,2)+oldmesh.nbs];
newmesh.elm_som(4*(0:oldmesh.nbt-1)+3,:) = [oldmesh.elm_fac(1:oldmesh.nbt,2)+oldmesh.nbs ...
    oldmesh.elm_som(1:oldmesh.nbt,1) ...
    oldmesh.elm_fac(1:oldmesh.nbt,3)+oldmesh.nbs];
newmesh.elm_som(4*(0:oldmesh.nbt-1)+4,:) = [oldmesh.elm_fac(1:oldmesh.nbt,3)+oldmesh.nbs ...
    oldmesh.elm_som(1:oldmesh.nbt,2) ...
    oldmesh.elm_fac(1:oldmesh.nbt,1)+oldmesh.nbs];

% Coord et numero de zones des anciens sommets
newmesh.som_coo(1:oldmesh.nbs,:) = oldmesh.som_coo(1:oldmesh.nbs,:);
newmesh.som_zon(1:oldmesh.nbs  ) = oldmesh.som_zon(1:oldmesh.nbs  );

% Coord et numero de zones des nouveaux sommets
for nf=1:oldmesh.nba
    newmesh.som_coo(nf+oldmesh.nbs,1) = ...
        sum(oldmesh.som_coo(oldmesh.fac_som(nf,:),1))/2.0;
    newmesh.som_coo(nf+oldmesh.nbs,2) = ...
        sum(oldmesh.som_coo(oldmesh.fac_som(nf,:),2))/2.0;
    newmesh.som_zon(nf+oldmesh.nbs) = ...
        min(oldmesh.som_zon(oldmesh.fac_som(nf,:)));
end
        
% On se preoccupe maintenant des aretes de bords
% calcul des aretes de bord
next = [2 3 1];

newmesh.nbab = 0;
for ie = 1:newmesh.nbt
    for i = 1:nb_som_elm
        is1=newmesh.elm_som(ie,i);
        is2=newmesh.elm_som(ie,next(i));
        iz1=newmesh.som_zon(is1);
        iz2=newmesh.som_zon(is2);
        if (iz1~=0) & (iz2~=0) 
            newmesh.nbab = newmesh.nbab + 1;
            newmesh.abd_som(newmesh.nbab,1:2) = [is1 is2];
        end
    end
end

%-------------------------------------------------------------
% Programmes de Jan Valdman
% http://sites.google.com/site/janvaldman/

function [element2edges, edge2nodes, edge2elements]=getEdges(elements)
%function: [element2edges, edge2nodes]=edge_numbering(elements)
%requires: deleterepeatedrows
%generates edges of (triangular) triangulation defined in elements
%elements is matrix, whose rows contain numbers of its element nodes 
%element2edges returns edges numbers of each triangular element
%edge2nodes returns two node numbers of each edge
%example: [element2edges, edge2nodes]=edge_numbering([1 2 3; 2 4 3])

%extracts sets of edges 
edges1=elements(:,[2 3]);
edges2=elements(:,[3 1]);
edges3=elements(:,[1 2]);

%as sets of their nodes (vertices)
vertices=zeros(size(elements,1)*3,2);
vertices(1:3:end,:)=edges1;
vertices(2:3:end,:)=edges2;
vertices(3:3:end,:)=edges3;

%repeated sets of nodes (joint edges) are eliminated 
[edge2nodes,element2edges]=deleterepeatedrows(vertices);
element2edges=reshape(element2edges,size(elements,2),size(elements,1))';

edge2elements=entryInWhichRows(element2edges);

%-------------------------------------------------------------
function [matrix,I]=deleterepeatedrows(matrix)
%function: [element2edges, edge2nodes]=edge_numbering(elements)
%requires: deleterepeatedrows
%generates edges of (triangular) triangulation defined in elements
%elements is matrix, whose rows contain numbers of its element nodes 
%element2edges returns edges numbers of each triangular element
%edge2nodes returns two node numbers of each edge
%example: [element2edges, edge2nodes]=edge_numbering([1 2 3; 2 4 3])

[matrixs,tags] = sortrows(sort(matrix,2));

% which ones were reps?
k = find(all(diff(matrixs)==0,2));

%these rows of matrix are repeated 
repeated=tags(k);

%and these rows will be removed
removed=tags(k+1);

%both lists are sorted 
[removeds, tags2]=sort(removed);
repeateds=repeated(tags2);

% delete the tags to removed rows
tags(k+1) = [];
% and recover the original array, in the original order.
matrix = matrix(sort(tags),:);

%row indices before matrix compression indicating repetition   
I=insertvector((1:size(matrix,1))',repeateds,removeds);

%-------------------------------------------------------------
function r=insertvector(v,pos_from,pos_to)
     tf=false(1,numel(v)+numel(pos_to));
     r=double(tf);
     tf(pos_to)=true;
     r(~tf)=v;
     r(tf)=r(pos_from);    
     
%-------------------------------------------------------------
function entryrows=entryInWhichRows(A)
%function: entryrows=entryInWhichRows(A)
%requires: none
%for every entry of integer matrix A, 
%its rows indices are stored in output matrix,
%zeros entries indicate no more occurence
%example: entryrows=entryInWhichRows([1 2; 1 3; 2 2]) returns
%         entryrows=[1   2   0; 
%                    1   3   3; 
%                    2   0   0]   
%meaning: entry 1 appears in rows 1 and 2
%         entry 2 appears in rows 1 and 3 (twice)
%         entry 3 appears in row  2 only

%size computation; 
r=max(max(A));
repetition=accumarray(A(:),ones(numel(A),1));
c=max(repetition);

%filling rows occurences
%this part should be somehow vectorized!
entryrows=zeros(r,c);
repetition=zeros(r,1);
for i=1:size(A,1)
    for j=1:size(A,2)
       index=A(i,j); 
       repetition(index)=repetition(index)+1; 
       entryrows(index,repetition(index))=i;
    end
end
