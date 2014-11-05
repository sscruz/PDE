function problema()
% Resolucion del problema -u'' = f(x) (en x \in (extremoIzd, extremodrch))
% u(0) = alpha
% u(1) = beta
%
% uExacta: solucion exacta al problema
% N: numero de puntos interiores en la discretizacion 

%%%%%%%%%%%%%%%%
%%%% Datos: %%%%
%%%%%%%%%%%%%%%%
	extremoIzd = 0;
	extremoDrch = 1;
	
	alpha = 0;
	beta = 0;

	uExacta = @(x) x.*sin(pi*x);
	f = @(x) -2*pi*cos(pi*x) + pi*pi*x.*sin(pi*x);

	N = 200;

	h = (extremoDrch - extremoIzd)/(N+1);
	x = extremoIzd:h:extremoDrch;
	x = x';

	%puntos interiores
	xi = x(2:N+1);		


%%%%%%%%%%%%%%%%%%%
%%%% Programa: %%%%
%%%%%%%%%%%%%%%%%%%

	%Calculo de matrices y vectores
	A = tridiag(N)/(h*h);
	b = f(xi);
	b(1) = b(1) + alpha/(h*h);
	b(N) = b(N) + beta/(h*h);


	% resolucion del problema
	uAproximada = A\b;
	
	% calculo del error

%%%%%%%%%%%%%%%%%%
%%%% Dibujos: %%%%
%%%%%%%%%%%%%%%%%%
	% Limpiar dibujos previos
	close 

	%Dibujar
	plot(x,[alpha;uAproximada;beta],'*');
	axis([extremoIzd, extremoDrch, min(uAproximada)-0.2, max(uAproximada)+0.2])	
	hold on

	plot(x,uExacta(x),'r');

	hold off;
