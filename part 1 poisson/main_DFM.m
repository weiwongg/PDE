%---------PART 3: solve PDE when f=2*pi^2*sin(pi*x).*sin(pi*y) and compare error among different step-----------
f = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y);

g = @(x,y) sin(pi*x).*sin(pi*y); %analytical solution of this problem
step=-1:-1:-10;
h=2.^step;
n=1./h;
err=zeros(1,10);
for i=1:10;
	nn=n(i);
	x=(1:1:nn-1)./nn;
	y=(nn-1:-1:1)./nn;
	[X,Y] = meshgrid(x,y);
	F=h(i)^2*f(X,Y);
	F=F(:);
	A = fdm2D(nn);
	u_h=A\F; %estimated solution of this problem
	u=g(X,Y);
	u=u(:);
	err(1,i)=norm(u-u_h,'inf');
end
figure(1)
plot(err)%plot the error with respect to i
xlabel('i')
ylabel('error')
title('error decrease with respect to i')
