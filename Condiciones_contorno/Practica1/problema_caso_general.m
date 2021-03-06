% Resolucion del problema -a(x)*u''+b(x)*u'(x)+c(x)*u(x) = f(x) (en x \in (extremoIzd, extremodrch))
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
	
	alpha = 1;
	beta = 1;

	uExacta = @(x) exp(-5*x.*(1-x));
	f = @(x) 0*x;
	a = @(x) 1+0*x;
	b = @(x) 5*(2*x-1);
	c = @(x) 10 + 0*x;

	N = 100;

	h = (extremoDrch - extremoIzd)/(N+1);
	x = extremoIzd:h:extremoDrch;
	x = x';

	%puntos interiores
	xi = x(2:N+1);		


%%%%%%%%%%%%%%%%%%%
%%%% Programa: %%%%
%%%%%%%%%%%%%%%%%%%

	%Calculo de matrices y vectores
	A =  sparse(1:N,1:N,2*a(x(2:end-1))/(h*h)+c(x(2:end-1)),N,N);
	A += sparse(1:N-1,2:N,-a(x(3:end-1))/(h*h)-b(x(3:end-1))/(2*h),N,N)';
        A += sparse(1:N-1,2:N,-a(x(2:end-2))/(h*h)+b(x(2:end-2))/(2*h),N,N);

        R = f(xi);
        R(1) = R(1) + a(xi(1))*alpha/(h*h) + b(xi(1))*alpha/(2*h);
        R(N) = R(N) + a(xi(N))*beta/(h*h) - b(xi(N))*beta/(2*h);


	% resolucion del problema
	uAproximada = A\R;
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
