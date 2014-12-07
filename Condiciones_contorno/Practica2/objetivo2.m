% Resolucion del problema a(x)*u''+c(x)*u(x) = f(x) (en x \in (extremoIzd, extremodrch))
% u'(0) = alpha
% u(1) = beta
%
% uExacta: solucion exacta al problema
% N1: numero de puntos interiores en la discretizacion en (0,1/2)
% N2: numero de puntos interiores en la discretizacion en (1/2,1)
% 
%


%%%%%%%%%%%%%%%%
%%%% Datos: %%%%
%%%%%%%%%%%%%%%%
	extremoIzd = 0;
	extremoDrch = 1;
	puntoMedio = 1/2;
	
	alpha = 0;
	beta = 0;
	N1 = 10;
	N2 = 20;

	uExacta = @(x) (1-x).*exp(x);
	f = @(x) (1+2*x-x.*x).*exp(x);
	a = @(x) 1+0*x;
	c = @(x) x;



	
%%%%%%%%%%%%%%%%%%%
%%%% Programa: %%%%
%%%%%%%%%%%%%%%%%%%
	%Discretizacion
	N = N1+N2+1;
	h1 = (puntoMedio - extremoIzd)/(N1+1);
	h2 = (extremoDrch - puntoMedio)/(N2+1);
	
	x1 = extremoIzd:h1:puntoMedio;
	x2 = puntoMedio:h2:extremoDrch;
	
	x = [x1(1:end-1),x2]';
	xi = x(2:end-1);
	h = x(2:end)-x(1:end-1);
	h = [h(1);h]; %vector de pasos extendido para hacer la matriz
		       %correctamente, defino h0 = h1

	%Calculo de matrices
	A =  sparse(1:N+1,1:N+1,2*a(x(1:end-1))./(h(1:end-1).*h(2:end))+c(x(1:end-1)),N+1,N+1);
	A += sparse(1:N,2:N+1,-2*a(x(2:end-1))./(h(2:end-1).*(h(2:end-1)+h(3:end))),N+1,N+1)';	
	A += sparse(1:N,2:N+1,-2*a(x(1:end-2))./(h(2:end-1).*(h(2:end-1)+h(1:end-2))),N+1,N+1);
	A(1,1) = A(1,1)/2;
	
	F = f(xi);
	F = [f(extremoIzd)/2-alpha/h(1);F];
	
	uAproximada = A\F;

%%%%%%%%%%%%%%%%%%
%%%% Dibujos: %%%%
%%%%%%%%%%%%%%%%%%

	% Limpiar dibujos previos
	close 

	%Dibujar

	%Dibujar
	plot(x,[uAproximada;beta],'*');
	axis([extremoIzd, extremoDrch, min(uAproximada)-0.2, max(uAproximada)+0.2])	
	hold on

	plot(x,uExacta(x),'r');

	hold off;

