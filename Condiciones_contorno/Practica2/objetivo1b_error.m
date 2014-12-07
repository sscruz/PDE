% Resolucion del problema a(x)*u''+c(x)*u(x) = f(x) (en x \in (extremoIzd, extremodrch))
% u'(0) = alpha
% u(1) = beta
%
% uExacta: solucion exacta al problema
% N: numero de puntos interiores en la discretizacion 
%
% Se desarrolla el caso del objetivo 1
% Se utiliza la formula u'(t) = (u(t+h)-u(t)) / h para aproximar la
% condicion de contorno en x = 0,
%



%%%%%%%%%%%%%%%%
%%%% Datos: %%%%
%%%%%%%%%%%%%%%%
	extremoIzd = 0;
	extremoDrch = 1;
	
	alpha = 0;
	beta = 0;

	uExacta = @(x) (1-x).*exp(x);
	f = @(x) (1+2*x-x.*x).*exp(x);
	a = @(x) 1+0*x;
	c = @(x) x;

	Nvec = 2.^(4:20)-1; % distintos valores para N
	err = [];



%%%%%%%%%%%%%%%%%%%
%%%% Programa: %%%%
%%%%%%%%%%%%%%%%%%%
	for N = Nvec,
	h = (extremoDrch - extremoIzd)/(N+1);
	x = extremoIzd:h:extremoDrch;
	x = x';

	%puntos interiores
	xi = x(2:N+1);		

	%Calculo de matrices y vectores
 
	A =  sparse(1:N+1,1:N+1,2*a(x(1:end-1))/(h*h)+c(x(1:end-1)),N+1,N+1);
	A += sparse(1:N,2:N+1,-a(x(2:end-1))/(h*h),N+1,N+1)';
        A += sparse(1:N,2:N+1,-a(x(1:end-2))/(h*h),N+1,N+1);
	A(1,1) = A(1,1)/2;

        R = f(xi);
        
	R = [f(extremoIzd)/2-alpha*a(extremoIzd)/h;R];
        R(end) += a(xi(N))*beta/(h*h);


	% resolucion del problema
	uAproximada = A\R;

	% calculo del error
	err = [err,max(abs(uAproximada-uExacta(x(1:end-1))))];

	end
%%%%%%%%%%%%%%%%%%
%%%% Dibujos: %%%%
%%%%%%%%%%%%%%%%%%
	% Limpiar dibujos previos
	close 
	loglog((extremoDrch - extremoIzd)./(Nvec+1),err,'+r');
	hold on
	loglog((extremoDrch - extremoIzd)./(Nvec+1),((extremoDrch - extremoIzd)./(Nvec+1)).^2,'b')

	%Dibujar


