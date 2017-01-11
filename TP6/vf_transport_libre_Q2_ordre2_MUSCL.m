clear all; 

% Creation de la fonction minmod pour implementer MUSCL :
beta = 2.0;
minmod = @(a,b) (sign(a).*max(0,sign(a.*b).*min(abs(a),abs(b))));
sweby = @(a,b) (sign(a).*max(0,sign(a.*b)).*max( min(beta*abs(a),abs(b)) , min(abs(a),beta*abs(b)) ) );

% Donnees du probleme :
a = 1;
N = 1600;
h = 1.0/N;
nu = .1; % <1/2
dt = nu*h/a;
T=10;

% Definition des vecteurs utiles la vectorisation :
LEFT = [N, 1:N];
RIGHT = [1:N, 1];

% Vecteur x du centre des mailles : 
xj = h/2:h:1-h/2;

% Definition de la fonction qui permet de dériver le type de schéma
% (upwind, LFM, LW, De Vuyst-Jaisson) à partir du schéma générique
varphi = @(nu) (1); % upwind
%varphi = @(nu) (1.0/nu); % LFM
%varphi = @(nu) (nu); %LW
%varphi = @(nu) (sqrt(nu)); %De Vuyst-Jaisson (1)
%varphi = @(nu) (nu + .25*(1-(2*nu-1)^2)); % De Vuyst-Jaisson (2)

% Donnee initiale : 
u0 = max(sin(6*pi*xj),0).*(xj<=1/3) ...
   + (3*xj-1).*(xj>1/3).*(xj<2/3) ...
   + (1).*(xj >=2/3);

clf(); plot(xj,u0, 'o-'); grid;
drawnow;

% initialisation 
t = 0;
u = u0;
i = 0;

% Boucle sur les pas de temps
while(t < T)
    i = i+1;
    t = t+dt;
    
    
    % calcul des pentes
    %sj = minmod(u(RIGHT(1:N))-u,u-u(LEFT(1:N))); % minmod 
    sj = sweby(u([2:N,1])-u,u-u([N,1:N-1])); % sweby
    % calul des valeurs en bords d'intervalle
    uplus = u + .5*sj;
    umoins = u - .5*sj;
    
    % Flux numerique
    phi = .5*a*(uplus(LEFT)+umoins(RIGHT)) ...
        -.5*abs(a)*varphi(nu) ...
        * (umoins(RIGHT) - uplus(LEFT));
    u = u - dt/h*(phi(2:N+1) - phi(1:N));
    if mod(i,10) == 0
       clf(); plot(xj,u, 'o-'); grid;
       drawnow;
    end
end
    

