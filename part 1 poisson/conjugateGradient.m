function [u, iter_num] = conjugateGradient(A, f, u0, tol)
	d=f-A*u0;
	r=d;
	u=u0;
	err=norm(r,2)^2;
	iter_num=1;
	while(sqrt(err)>tol)
		z=A*d;
		alpha=err^2/(d'*z);
		u=u+alpha*d;
		r = r-alpha*z;
		err_old=err;
		err=norm(r,2)^2;
		beta=err/err_old;
		d = r + beta*d;
		iter_num = iter_num + 1;
	end
end
