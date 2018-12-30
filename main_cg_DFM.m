%---------PART 4: using conjugate gradient descent method to solve PDE---------------
ff = @(x,y) exp(-(x-y).^2); 
hh=2^(-10); %step
nn=2^10;
x=(1:1:nn-1)/nn;
y=(nn-1:-1:1)./nn;
[X,Y] = meshgrid(x,y);
FF=hh^2*ff(X,Y);
FF=FF(:);
AA = fdm2D(nn);
uu_h=AA\FF; %solution using mldivide
u0=zeros((nn-1)^2,1);
[uu_hConj, iter_num] = conjugateGradient(AA, FF, u0,1e-8); %solution using conjugate gradient descent method

step=-1:-1:-10;
h=2.^step;
n=1./h;
iters_num=zeros(1,10);

for i=1:10;
	nn=n(i);
	x=(1:1:nn-1)./nn;
	y=(nn-1:-1:1)./nn;
	[X,Y] = meshgrid(x,y);
	F=ff(X,Y);
	F=F(:);
	A = fdm2D(nn);
	[u_hConj, iter_num]=conjugateGradient(A, F, u0,1e-6);
	iters_num(1,i)=iter_num;
end
figure(2)
plot(iters_num)%plot the error with respect to i
xlabel('i')
ylabel('iteration number')
title('the number of iterations versus i')
