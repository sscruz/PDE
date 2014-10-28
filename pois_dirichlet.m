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


%%%%%%%%%%%%
%% Dibujo %%
%%%%%%%%%%%%


subplot(1,2,1);
surf(X,Y,uExacta(X,Y));
axis([a a+L b b+L -1 1]);

subplot(1,2,2);
surf(X,Y,reshape(uAprox,N+2,N+2));
axis([a a+L b b+L -1 1]);
