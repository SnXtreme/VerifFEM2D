function h_adapt = adapt(omega,error,error_0)
h=omega.elementSize;
p=omega.order;
d=omega.dim;
h_adapt=(h*error_0^(1/p))./(error.^(1/(2*p+d))*((sum(error.^(d/(2*p+d)))).^(1/(2*p))));
    