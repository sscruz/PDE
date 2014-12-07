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
	epsilon = .01;
	
	uExacta = @(x) (exp(x/epsilon)-1)./(exp(1/epsilon)-1);
	f = @(x) 0*x;
	b = 1;



	
%%%%%%%%%%%%%%%%%%%
%%%% Programa: %%%%
%%%%%%%%%%%%%%%%%%%
	%Discretizacion arbitraria
	x = linspace(0,.95,12);
	x = [x(1:end-1),linspace(0.95,1,200)];
	x = x';
	xi = x(2:end-1);
	N = length(xi);
	h = x(2:end)-x(1:end-1);
	
	
	%Calculo de matrices
	A = sparse(1:N,1:N,(h(1:end-1).*h(2:end)).\(2*epsilon)+h(1:end-1).\b,N,N);
        A +=sparse(1:N-1,2:N,-(h(2:end-1).*(h(2:end-1)+h(3:end))).\2*epsilon-h(2:end-1).\b,N,N)';
											
	A +=sparse(1:N-1,2:N,-(h(2:end-1).*(h(2:end-1)+h(1:end-2))).\2*epsilon,N,N);
	
	F = f(xi);
	F(1) += alpha*(2*epsilon/(h(1)*(h(1)+h(2)))+b/h(1));
	F(end) += (2*epsilon/(h(end)*(h(end)+h(end-1))))*beta;
	

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

