function [M]=face_number(M)
% M=face_number(M) modifie la structure maillage M 
%  
% Ce pgm (i)   numerote les aretes de 1 a M.nba 
%        (ii)  associe a chaque aretes les numeros des 2 elements 
%              adjacents (par ordre croissant des numeros) 
%              M.fac_elm(nf,i) == numero du ieme element adjacent 
%                                    a l'arete nf
%              NB : M.fac_elm(nf,2) == 0 si nf est une arete 
%                                           sur la frontiere du domaine
%       (iii)  associe a chaque arete nf les 2 numeros des sommets de
%              l'arete 
%              M.fac_som(nf,i) == numero du ieme sommet de l'arete 
%              Ces sommets sont ranges de telle sorte que la normale a 
%              la face orientee de M.fac_elm(nf,1) vers M.fac_elm(nf,2) 
%              soit egale a un coefficient > 0 pres a
%              [M.som_coo(M.fac_som(:,2),2) - M.som_coo(M.fac_som(:,1),2) ; 
%              -M.som_coo(M.fac_som(:,2),1) + M.som_coo(M.fac_som(:,1),1)].
%       (iv)   calcule un numero de zone pour chaque arete :
%              M.fac_zon(nf)  
%       (v)    calcule pour chaque arete le milieu de l'arete :
%              M.fac_gra(nf,1:2)
%       (vi)   calcule pour chaque arete nf la normale unitaire de l'arete orientee 
%              de M.fac_elm(nf,1) vers M.fac_elm(nf,2) :
%              M.fac_nor(nf,1:2)
%              M.fac_zon(nf)
%       (vii)  calcule la longueur de l'arete : 
%              M.fac_mes(nf)
%       (viii) calcule la surface de chaque triangle ie : 
%              M.elm_mes(ie) 
%       (ix)   calcule le centre de gravite de chaque triangle ie : 
%              M.elm_gra(ie) 
%        (x)   calcule les aretes de chaque element
%              M.elm_fac(ie,i) == numero de la ieme face de l'element ie
%
%Copyright (c) 2005 by Frederic Pascal, ENS de CACHAN
%Copyright (c) 2009 by Jan Valdam, http://sites.google.com/site/janvaldman/software 

% Verification
if isfield(M,'nbs') ~= 1 ...
        |  isfield(M,'nbt') ~= 1 ...
        |  isfield(M,'som_coo') ~= 1 ...
        |  isfield(M,'som_zon') ~= 1 ...
        |  isfield(M,'elm_som') ~= 1
    error([' LA STRUCTURE OU LES CHAMPS DE LA STRUCTURE M NE SONT PAS CORRECT :' ...
            ' RENOMMEZ VOS CHAMPS '])
end

% Verification
nb_som_elm = size(M.elm_som,2);
if ( nb_som_elm ~= 3)
    error([' NB DE SOMMETS PAR VOLUME DIFFERENT DE 3 '])
end 

% Appel aux programmes de Jan Valdam 
[M.elm_fac,M.fac_som,M.fac_elm]=getEdges(M.elm_som);

% Zone des faces
M.fac_zon = min(M.som_zon(M.fac_som(:,:)),[],2);

% Nba 
M.nba = size(M.fac_som,1);

% Calcul du centre de gravite des faces et normales * mesure des faces
coo_isf(:,:,1)=M.som_coo(M.fac_som(:,1),:);
coo_isf(:,:,2)=M.som_coo(M.fac_som(:,2),:);

M.fac_gra= 0.5*sum(coo_isf,3);

veca(:,:)=coo_isf(:,:,2)-coo_isf(:,:,1);

M.fac_mes=sqrt(veca(:,2).^2+veca(:,1).^2);
M.fac_nor=[veca(:,2)./M.fac_mes,-veca(:,1)./M.fac_mes];

clear veca
clear coo_isf

% calcul de la surface de ie et du centre de gravite de ie
coo_ise(:,:,1) = M.som_coo(M.elm_som(:,1),:);
coo_ise(:,:,2) = M.som_coo(M.elm_som(:,2),:);
coo_ise(:,:,3) = M.som_coo(M.elm_som(:,3),:);

M.elm_gra=sum(coo_ise,3)/3;

vecb(:,:,1:2) = coo_ise(:,:,[1 2]) - coo_ise(:,:,[3 3]);
M.elm_mes = 0.5*(vecb(:,1,1).*vecb(:,2,2)-vecb(:,1,2).*vecb(:,2,1));

clear vecb

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
