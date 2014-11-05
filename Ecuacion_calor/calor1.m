%% Resolucion numerica de la eq del calor
%% u_t (x,t) -nu*u_xx(x,t) = f(x,t)
%% u(


%%%%%%%%%%%
%% Datos %%
%%%%%%%%%%%
    
a = 0;
b = 1;
T =2;

%uInicial = @(x) 4*(x-a).*(b-x);
uInicial = @(x) 10*(0.25-x).*(x-0.75).*(x>0.25).*(x<0.75);
alpha = @(t) sin(4*pi*t);
beta = @(t) 4*t.*(T-t);
nu = 1;
fuente = @(x,t) x.*t;



N = 63;
M = 50;

%%%%%%%%%%%%%%
%% Programa %%
%%%%%%%%%%%%%%
h = (b-a)/(N+1);
tau = T/M;
x = [a:h:b]';
xi = x(2:N+1);




lambda = tau*nu/(h*h);

D = tridiag(N);
A = speye(N) + lambda*D;

u0 = uInicial(xi);
AllU = u0;
PlotU = [alpha(0);u0;beta(0)];

plot(x,PlotU);

for k=1:M,
b = tau*fuente(xi,tau*k);
b(1) = lambda*alpha(k*tau);
b(end) = lambda*beta(k*tau);

b += AllU(:,k);
AllU(:,k+1) = A\b;
PlotTmp = [alpha(k*tau);AllU(:,k+1);beta(k*tau)];
PlotU = [PlotU, PlotTmp];

hold on
plot(x,PlotU(:,k+1));
pause(.1);

end




%%%%%%%%%%%%
%% Dibujo %%
%%%%%%%%%%%%

tmp = 0:tau:T';
[X,TMP] = meshgrid(x,tmp);
hold off
close

waterfall(X',TMP',PlotU);
