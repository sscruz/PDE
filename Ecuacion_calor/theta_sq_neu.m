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

theta = 1; %coeficiente del theta-esquema

%uInicial = @(x) 4*(x-a).*(b-x);
%uInicial = @(x) 10*(0.25-x).*(x-0.75).*(x>0.25).*(x<0.75);


%alpha = @(t) sin(4*pi*t); %temperatura en x = a
%beta = @(t) 4*t.*(T-t); %temperatura en x = b
alpha = @(t) 0*t;
beta = @(t) 0*t;

nu = 1;
fuente = @(x,t) 0*x+0*t;



N = 100; % subdivisiones en el espacio
M = 40; % subdivisiones en el tiempo

%%%%%%%%%%%%%%
%% Programa %%
%%%%%%%%%%%%%%
h = (b-a)/(N+1);
tau = T/M;
x = [a:h:b]';


lambda = tau*nu/(h*h);

D = tridiag(N+2);
D(1,2) = -2.;
D(end,end-1) = -2.;

A = speye(N+2) + theta*lambda*D;
[V,U] = lu(A); % descomposicion LU de la matriz A

B = speye(N+2) - (1-theta)*lambda*D;

u0 = uInicial(x);
AllU = 20*(abs(x-(b+a)/2) <= (b-a)/4);

plot(x,AllU);

b = tau*fuente(x,0);
b(1) = lambda*alpha(0);
b(end) = lambda*beta(0);
  
%Resolucion para cada valor del tiempo
for k=1:M,
    b = [b,tau*fuente(x,tau*k)];
%Actualizar esto para las cond. von Neumann 
    b(1,k+1) = 0; %lambda*alpha(k*tau);
    b(end,k+1) = 0;% lambda*beta(k*tau);
    
    c = theta*b(:,k+1)+(1-theta)*b(:,k) + B*AllU(:,k);
    
    AllU(:,k+1) = A\c;
  
    hold on
    plot(x,AllU(:,k+1));
    pause(.1);
  
end




%%%%%%%%%%%%
%% Dibujo %%
%%%%%%%%%%%%

tmp = 0:tau:T';
[X,TMP] = meshgrid(x,tmp);
hold off
close

waterfall(X',TMP',AllU);
