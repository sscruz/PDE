function error2()
% Calculo del error del problema:
% -u'' = f(x) (en x \in (extremoIzd, extremodrch))
% u'(0) = alpha
% u(1) = beta
%
% uExacta: solucion exacta al problema
% Nvec: vector de numero de puntos en la discretizacion para los que
% se calcula el error

%%%%%%%%%%%%%%%%
%%%% Datos: %%%%
%%%%%%%%%%%%%%%%
	extremoIzd = 0;
	extremoDrch = 1;
	
	alpha = 1;
	beta = -1;

	uExacta = @(x) x.*cos(pi*x);
	f = @(x) 2*pi*sin(pi*x) + pi*pi*x.*cos(pi*x);

	Nvec = 2.^(4:15)-1; % distintos valores para N
	error = [];

%%%%%%%%%%%%%%%%%%%
%%%% Programa: %%%%
%%%%%%%%%%%%%%%%%%%

	for N = Nvec

		h = (extremoDrch - extremoIzd)/(N+1);
		x = extremoIzd:h:extremoDrch;
		x = x';
	
		%puntos interiores que se resuelven
		xi = x(1:N+1);		



		%Calculo de matrices y vectores
		A = tridiag(N+1);
		A(1,1) = 1;
		A = A/(h*h);
		b = f(xi);
		b(1) = -alpha/h;
		b(N+1) = b(N+1) + beta/(h*h);
	

		% resolucion del problema
		uAproximada = A\b;
		
		% calculo del error
		error = [error; max(abs(uAproximada-uExacta(xi)))];	
	end
	
%%%%%%%%%%%%%%%%%%
%%%% Dibujos: %%%%
%%%%%%%%%%%%%%%%%%
	% Limpiar dibujos previos
	close 
	
	%Dibujar
	loglog((extremoDrch - extremoIzd)./(Nvec+1),error,'+r');
	hold on
	loglog((extremoDrch - extremoIzd)./(Nvec+1),((extremoDrch - extremoIzd)./(Nvec+1)).^2,'b')
