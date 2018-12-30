%----------PART 1: implementaion of F the discretised right hand side----------
function F = rhs(@(x,y) f(x,y), X, Y)
    F=f(X,Y);
	F=F(:);
end