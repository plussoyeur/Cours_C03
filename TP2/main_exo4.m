%  Resolution du systeme element fini de l'equation de Poisson avec
%  conditions aux limites de Dirichlet et de Neumann : 
%         - Du = f = -4 sur [0,1]*[0,1]
%         u = g1 = x² + y² sur GammaD
%         grad(u).n = 0 sur GammaN
%         
clear all;
mesh = lect_mesh('car0');




for k = 1:4
    
    h(k) = min(mesh.som_coo((find(mesh.som_coo(:,1) ~= 0 & mesh.som_zon ~= 0))));
    
    % Coefficient de diffusion
    kappa = ones(mesh.nbt,1);
    
    % Declaration de f
    f = @(x,y) (-4); % fonction anonyme
    
    % Declaration de g
    g = @(z,x,y) ((z==2).*(x.^2+y.^2));
    
    % Assemblage des matrices
    A = assemb_A(kappa,mesh);
    F = assemb_F(f,mesh);
    M = assemb_M(mesh);
    
    % Initialisation de u
    u = zeros(mesh.nbs,1);
    
    % Assemblage du second membre complet
    x = mesh.som_coo(:,1);
    y = mesh.som_coo(:,2);
    dir = find(mesh.som_zon == 2);
    inconnues = setdiff((1:mesh.nbs)',dir);
    u(dir) = g(mesh.som_zon(dir),x(dir),y(dir));
    F = F - A*u;
    
    % Resolution du systeme
    u(inconnues) = A(inconnues, inconnues)\F(inconnues);
    
    % Representation graphique
    %tri = mesh.elm_som;
    %figure(1); clf; colormap(cool);
    %trimesh(tri,x,y,u);
    
    % Vraie solution du problème
    u_real = x.^2 + y.^2;
    diff = u_real - u;
    
    % Erreur
    err(k) = sqrt(diff'*M*diff);
    
    if (k~=4)
        mesh = raf_mesh(mesh);
    end
    
    
end


figure();
loglog(h, err, 'rx');
hold on;
t = min(h):(max(h)-min(h))/100:max(h);
loglog(t,t*err(1)/h(1),'-','Color','blue');
loglog(t,t.^2*err(1)/h(1)^2,'-','Color','green');
loglog(t,t.^3*err(1)/h(1)^3,'-','Color','magenta');
loglog(t,t.^4*err(1)/h(1)^4,'-','Color','cyan');


