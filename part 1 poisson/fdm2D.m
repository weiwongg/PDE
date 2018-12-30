%----------PART 2: implementation of finite difference system matrix in two spatial dimensions---------
function A = fdm2D(n)
	e1=ones(n-1,1);
    e2=ones(n-2,1);
	T=diag(2*e1,0)+diag(-e2,-1)+diag(-e2,1);
    T=sparse(T);
	I=speye(n-1);
	A=kron(I,T)+kron(T,I);
end