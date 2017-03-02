function [ u0 ] = init(mesh)
%Emeriau Pierre-Emmanuel
%emeriau.pe@gmail.co

f = @(x,y) max( 4*(0.25 - sqrt( (x-.5).^2 + y.^2)  ), 0);
x = mesh.elm_gra(:,1);
y = mesh.elm_gra(:,2);

u0 = zeros(mesh.nbt);
u0 = f(x,y);


end

