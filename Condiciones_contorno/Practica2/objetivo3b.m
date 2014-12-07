% Resolucion del problema -epsilon*u''+b*u'(x) = f(x) (en x \in (extremoIzd, extremodrch))
% u(extremoIzd) = alpha
% u(extremoDrch) = beta
% epsilon, b constantes
%
% uExacta: solucion exacta al problema
% 
%


%%%%%%%%%%%%%%%%
%%%% Datos: %%%%
%%%%%%%%%%%%%%%%

	extremoIzd = 0;
	extremoDrch = 1;
	
	alpha = 0;
	beta = 1;
	N = 10;
	epsilon = .01;
	
	uExacta = @(x) (exp(x/epsilon)-1)./(exp(1/epsilon)-1);
	f = @(x) 0*x;
	b = 1;



	
%%%%%%%%%%%%%%%%%%%
%%%% Programa: %%%%
%%%%%%%%%%%%%%%%%%%
	%Discretizacion
	h = (extremoDrch-extremoIzd)/(N+1);
	x = extremoIzd:h:extremoDrch;
	x = x';
	xi = x(2:end-1);

	
	%Calculo de matrices
	A = epsilon*tridiag(N)/(h*h);
	A += sparse(1:N,1:N,b/h,N,N);
	A += sparse(1:N-1,2:N,-b/h,N,N)';
	
	F = f(xi);
	F(1) += (epsilon/(h*h)+b/h)*alpha;
	F(end) += (epsilon/(h*h))*beta;

	uAproximada = A\F;

%%%%%%%%%%%%%%%%%%
%%%% Dibujos: %%%%
%%%%%%%%%%%%%%%%%%

	% Limpiar dibujos previos
	close 

	%Dibujar

	%Dibujar
	plot(x,[alpha;uAproximada;beta],'*');
	axis([extremoIzd, extremoDrch, min(uAproximada)-0.2, max(uAproximada)+0.2])	
	hold on

	plot(x,uExacta(x),'r');

	hold off;

