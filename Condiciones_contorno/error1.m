function error1()
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

	Nvec = 2.^(4:20)-1; % distintos valores para N
	error = [];

%%%%%%%%%%%%%%%%%%%
%%%% Programa: %%%%
%%%%%%%%%%%%%%%%%%%

	for N = Nvec

		h = (extremoDrch - extremoIzd)/(N+1);
		x = extremoIzd:h:extremoDrch;
		x = x';
	
		%puntos interiores
		xi = x(2:N+1);		



		%Calculo de matrices y vectores
		A = tridiag(N)/(h*h);
		b = f(xi);
		b(1) = b(1) + alpha/(h*h);
		b(N) = b(N) + beta/(h*h);
	

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
