function x = newton_ralph(func,smallstep,ytol,xguess,Rootval)
x=xguess;
yerr=Rootval-func(x);
while abs(yerr)>ytol
    dydx=(func(x+smallstep)-func(x))/smallstep;
    dx=yerr/dydx;
    x=x+dx;
    yerr = Rootval-func(x);
end