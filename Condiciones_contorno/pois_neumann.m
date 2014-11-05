%%%%%%%%%%%
% Resolucion del problema u_xx + u_yy = f(x,y) en un cuadrado con vertices 
% (a,b) -- (a+L,b) -- (a+L,b+L) -- (a,b+L)
% y condiciones de contorno mixtas
% N: tamano de la discretizacion


% Nota: se puede variar que condiciones son de tipo Dirichlet y cuales
% Neumann comentando las lineas 47 a 61
% qO, qN, qS, qE son las condiciones tipo Neumann para las fronteras
% Oeste, Norte, Sur y Este.

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
qO = @(x,y) pi*sin(pi*x).*cos(pi*y); 
qN = @(x,y) -pi*cos(pi*x).*sin(pi*y);
qS = @(x,y) pi*cos(pi*x).*sin(pi*y);
qE = @(x,y) -pi*sin(pi*x).*cos(pi*y); 


N = 39;
h = L/(N+1);

%%%%%%%%%%%
%% Malla %%
%%%%%%%%%%%

x = a:h:a+L;
y = b:h:b+L;

[X,Y] = meshgrid(x,y);

xx = X(:);
yy = Y(:);


n_dirich = find ( abs(xx - a) < eps...
    | abs(yy - (b + L)) < eps);

n_neuS = find( abs(yy - b) < eps); % sur
n_neuE = find( abs(xx - (a + L)) < eps); % este


%%%%%%%%%%%%%%
%% Programa %%
%%%%%%%%%%%%%%

R = f(xx,yy);
R(n_dirich) = uExacta(xx(n_dirich),yy(n_dirich))/eps;
R(n_neuE) += 2*qE(xx(n_neuE),yy(n_neuE))/h;
R(n_neuS) += 2*qS(xx(n_neuS),yy(n_neuS))/h;


D = tridiag(N+2)/(h*h);
D(1,2) = -2/(h*h);
D(N+2,N+1) = -2/(h*h);

A = kron(D,speye(N+2)) + kron(speye(N+2),D);

for i = 1:length(n_dirich),
   A(n_dirich(i),n_dirich(i)) = 1/eps;
end

uAprox = A\R;


%%%%%%%%%%%%
%% Dibujo %%
%%%%%%%%%%%%


subplot(1,2,1);
surf(X,Y,uExacta(X,Y));
axis([a a+L b b+L -1 1]);

subplot(1,2,2);
surf(X,Y,reshape(uAprox,N+2,N+2));
axis([a a+L b b+L -1 1]);
