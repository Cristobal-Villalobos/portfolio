function nx=norm21(x)
%here we are calculating the l1l2 norm of x a magnetization ordered:
%x1
%x2
%x3.

    n=size(x,1)/3;

    x1  =x(1:n);
    x2  =x(n+1:2*n);
    x3  =x(2*n+1:3*n);
    xm  =[x1,x2,x3];
    nx  =sum(sqrt(sum(abs(xm).^2,2)));
    
    
end