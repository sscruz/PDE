function A=tridiag(N)
% Construye una matriz tridiagonal.
% 	2 en la diagonal principal
%	-1 en las primeras diagonales
% Argumentos:
%	+ N: dimension de la matriz

	E = sparse(1:N-1,2:N,ones(N-1,1),N,N);
	A = 2*speye(N,N) - E - E';

	
