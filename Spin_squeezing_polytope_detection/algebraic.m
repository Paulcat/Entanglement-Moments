function [h1] = algebraic(N)


%First order momenta is fixed:

x0=0;
y0=0;
z0=0;


%Some constant
k=(N-1)/N;


%Compute the vertices
Ax=[N^2/4-k.*(y0^2+z0^2), N/4+k.*(y0)^2, N/4+k.*(z0)^2];
Bx=[x0^2+(y0^2+z0^2)/N, N/4+k.*(y0)^2, N/4+k.*(z0)^2];

Ay=[N/4+k.*(x0)^2, N^2/4-k.*(x0^2+z0^2), N/4+k.*(z0)^2];
By=[N/4+k.*(x0)^2, y0^2+(x0^2+z0^2)/N, N/4+k.*(z0)^2];

Az=[N/4+k.*(x0)^2, N/4+k.*(y0)^2, N^2/4-k.*(x0^2+y0^2)];
Bz=[N/4+k.*(x0)^2, N/4+k.*(y0)^2, z0^2+(x0^2+y0^2)/N];

vertexs=[Ax; Ay; Az; Bx; By; Bz];

%Now do the polytope of the spin squeezing inequalities via alphaShape
h1 = alphaShape(vertexs, inf);


end