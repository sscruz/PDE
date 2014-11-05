%%%%%%%%%
%%Datos%%
%%%%%%%%%

%extremos 
a = 0;
b = 1;
L = 2;

%dato y soluciones 
uExacta = @(x,y) cos(pi*x).*cos(pi*y);
f = @(x,y) 2*pi*pi*cos(pi*x).*cos(pi*y);
qO = @(x,y) pi*sin(pi*x).*cos(pi*y); 
qN = @(x,y) -pi*cos(pi*x).*sin(pi*y);
qS = @(x,y) pi*cos(pi*x).*sin(pi*y);
qE = @(x,y) -pi*sin(pi*x).*cos(pi*y); 

Nvec = 2.^(4:8);
hvec = (Nvec+1).\L;

err = [];

for N = Nvec,
	h = L/(N+1);
	%%%%%%%%%%%
	%% Malla %%
	%%%%%%%%%%%

	x = a:h:a+L;
	y = b:h:b+L;

	[X,Y] = meshgrid(x,y);

	xx = X(:);
	yy = Y(:);


	n_dirich = find ( abs(xx - a) < eps...
	    | abs(yy - (b + L)) < eps);

	n_neuS = find( abs(yy - b) < eps); % sur
	n_neuE = find( abs(xx - (a + L)) < eps); % este


	%%%%%%%%%%%%%%
	%% Programa %%
	%%%%%%%%%%%%%%

	R = f(xx,yy);	
	R(n_dirich) = uExacta(xx(n_dirich),yy(n_dirich))/eps;
	R(n_neuE) += 2*qE(xx(n_neuE),yy(n_neuE))/h;
	R(n_neuS) += 2*qS(xx(n_neuS),yy(n_neuS))/h;


	D = tridiag(N+2)/(h*h);
	D(1,2) = -2/(h*h);
	D(N+2,N+1) = -2/(h*h);

	A = kron(D,speye(N+2)) + kron(speye(N+2),D);

	for i = 1:length(n_dirich),
   		A(n_dirich(i),n_dirich(i)) = 1/eps;
	end

	uAprox = A\R;
	
	err = [err,max(abs(uAprox-uExacta(xx,yy)))];

end

%%%%%%%%%%%%
%% Dibujo %%
%%%%%%%%%%%%
	close
	loglog(hvec,err,'*r');
	hold on
	loglog(hvec,hvec.^2);
	hold off
