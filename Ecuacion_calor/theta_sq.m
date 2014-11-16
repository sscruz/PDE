%% Resolucion numerica de la eq del calor usando met. Euler implicito
%% u_t (x,t) -nu*u_xx(x,t) = fuente(x,t) en a<x<b t<T
%% u(a,t) = alpha(t)
%% u(b,t) = beta(t)
%% u(x,0) = uInicial(x)


%%%%%%%%%%%
%% Datos %%
%%%%%%%%%%%

a = 0;
b = 1;
T = 1;

theta = 0.5; %coeficiente del theta-esquema

uInicial = @(x) 4*(x-a).*(b-x);
%uInicial = @(x) 10*(0.25-x).*(x-0.75).*(x>0.25).*(x<0.75);

alpha = @(t) sin(4*pi*t); %temperatura en x = a
beta = @(t) 4*t.*(T-t); %temperatura en x = b
nu = 1;
fuente = @(x,t) x.*t;		



N = 100; % subdivisiones en el espacio
M = 100; % subdivisiones en el tiempo

%%%%%%%%%%%%%%
%% Programa %%
%%%%%%%%%%%%%%
h = (b-a)/(N+1);
tau = T/M;
x = [a:h:b]';
xi = x(2:N+1);




lambda = tau*nu/(h*h);

D = tridiag(N);
A = speye(N) + theta*lambda*D;
B = speye(N) - (1-theta)*lambda*D;

u0 = uInicial(xi);
AllU = u0;
PlotU = [alpha(0);u0;beta(0)];

plot(x,PlotU);

d = tau*fuente(xi,0);
d(1) += lambda*alpha(0);
d(end) += lambda*beta(0);
  
%Resolucion para cada valor del tiempo
for k=1:M,
    d = [d,tau*fuente(xi,tau*k)];
    d(1,k+1) += lambda*alpha(k*tau);
    d(end,k+1) += lambda*beta(k*tau);
  
    c = theta*d(:,k+1)+(1-theta)*d(:,k) + B*AllU(:,k);
    
    AllU(:,k+1) = A\c;

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
