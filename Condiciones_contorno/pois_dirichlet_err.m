%%%%%%%%%%%
% Calculo del error para el problema 
% u_xx + u_yy = f(x,y) en un cuadrado con vertices 
% (a,b) -- (a+L,b) -- (a+L,b+L) -- (a,b+L)
% y condiciones de contorno Dirichlet
% N: tamano de la discretizacion

%%%%%%%%%
%%Datos%%
%%%%%%%%%

%extremos 
a = 0;
b = 1;
L = 2;

%dato y soluciones 
uExacta = @(x,y) cos(pi*x).*cos(pi*y);
f = @(x,y) 2*pi*pi*cos(pi*x).*cos(pi*y);

Nvec = 2.^(4:8);
hvec = (Nvec+1).\L;

err = [];

for N = Nvec
    %%%%%%%%%%  
    %% Malla %
    %%%%%%%%%%
    h = L/(N+1);

    x = a:h:a+L;
    y = b:h:b+L;

    [X,Y] = meshgrid(x,y);

    xx = X(:);
    yy = Y(:);


    nf = find ( abs(xx - a) < eps...
    | abs(xx - (a + L)) < eps |...
     abs(yy - b) < eps |...
     abs(yy - (b + L)) < eps);


    %%%%%%%%%%%%%%
    %% Programa %%
    %%%%%%%%%%%%%%


R = f(xx,yy);
R(nf) = uExacta(xx(nf),yy(nf))/eps;

D = tridiag(N+2)/(h*h);
A = kron(D,speye(N+2)) + kron(speye(N+2),D);

for i = 1:length(nf),
   A(nf(i),nf(i)) = 1/eps;
end

uAprox = A\R;

err = [err,max(abs(uAprox-uExacta(xx,yy)))];

end
%%%%%%%%%%%%
%% Dibujo %%
%%%%%%%%%%%%
close
loglog(hvec,err,'*r');
hold on
loglog(hvec,hvec.^2);
hold off
