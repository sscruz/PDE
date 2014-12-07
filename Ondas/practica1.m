% Resolucion de la ecuacion de ondas
% u_tt = c*c u_xx 
% u(a,t) = 0
% u(b,t) = 0
% u(x,0) = uInicial
% u_t(x,0) = dInicial

% Se toma una discretizacion en el intervalo [a,b] con N puntos
% interiores. 
% Se toma una discretizacion de M+1 puntos en el tiempo.
%

%%%%%%%
%Datos%
%%%%%%%


a = 0; % extremo izquierdo
b = 1; % extremo derecho
c = 1; % velocidad de la onda
T = 2; % tiempo total

% condiciones iniciales
%L = b-a;
%uInicial = @(x) x.*(L-x); % posicion inicial
%dInicial = @(x) 0*x; % velocidad inicial

%uInicial = @(x) ((x-0.25).*(-x+0.75) > 0).*(x-0.25).*(0.75-x);
%dInicial = @(x) 0*x; % velocidad inicial

uInicial = @(x) 0.1*sin(2*pi*x)+0.1*sin(pi*x)+0.1*sin(3*pi*x)+0.1*sin(4*pi*x);


%discretizacion espacial y temporal
N = 50; % numero de puntos en la discretizacion espacial
M = 1000; % numero de puntos en la discretizacion temporal

h = (b-a)/(N+1);  % paso espacial
x = (a:h:b)'; % nodos en los que se calcula la solucion
xi = x(2:end-1);

tau = T/M; % paso temporal


%%%%%%%%%%%%
% Programa %
%%%%%%%%%%%%

lambda = c*tau/h;

% test de estabilidad
if (lambda > 1) 
   printf("Metodo inestable");
   return;
end 

% construccion de matrices
A = -lambda*lambda*tridiag(N);
A += 2*speye(N);

% calculo de u(x,0) y u(x,t_1)
uAproximada = uInicial(xi);
uAproximada = [uAproximada, 0.5*A*uInicial(xi)+tau*dInicial(xi)];
plot(x,[0;uAproximada(:,1);0]);
pause(0.01);
plot(x,[0;uAproximada(:,2);0]);
pause(0.01);

% resolucion de la ecuacion en los distintos valores del tiempo
for k = 3:M+1,
    uAproximada(:,k) = A*uAproximada(:,k-1) - uAproximada(:,k-2);
    
    plot(x,[0;uAproximada(:,k);0]);
    axis([a,b,-.25,.25]);    
    pause(0.01);
end

%%%%%%%%%%
% Dibujo %
%%%%%%%%%%

% dibujo de todas las soluciones
tmp = 0:tau:T;
[X,TMP] = meshgrid(x,tmp);

close
waterfall(X',TMP',[zeros(1,M+1);uAproximada;zeros(1,M+1)]);
