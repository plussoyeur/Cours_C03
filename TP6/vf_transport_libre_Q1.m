clear all; 

% Donnees du probleme :
a = 1;
N = 1600;
h = 1.0/N;
nu = .1;
dt = nu*h/a;
T=4;

% Definition des vecteurs utiles la vectorisation :
LEFT = [N, 1:N];
RIGHT = [1:N, 1];

% Vecteur x du centre des mailles : 
xj = h/2:h:1-h/2;

% Definition de la fonction qui permet de dériver le type de schéma
% (upwind, LFM, LW, De Vuyst-Jaisson) à partir du schéma générique
%varphi = @(nu) (1); % upwind
%varphi = @(nu) (1.0/nu); % LFM
%varphi = @(nu) (nu); %LW
%varphi = @(nu) (sqrt(nu)); %De Vuyst-Jaisson (1)
varphi = @(nu) (nu + .25*(1-(2*nu-1)^2)); % De Vuyst-Jaisson (2)

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
    % Flux numerique
    phi = .5*a*(u(LEFT)+u(RIGHT)) ...
        -.5*abs(a)*varphi(nu) ...
        * (u(RIGHT) - u(LEFT));
    u = u - dt/h*(phi(2:N+1) - phi(1:N));
    if mod(i,10) == 0
       clf(); plot(xj,u, 'o-'); grid;
       drawnow;
    end
end
    

