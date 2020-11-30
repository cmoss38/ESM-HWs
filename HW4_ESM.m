%% HW 4 
%% 1. ...Produce a short Matlab script to test your answer.
sol=inv([7,2,1;0,3,-1;-3,4,-2]);

%% 2.a. Implement the bisection scheme to solve for E.
for i=1:1
w=.017; %set up values
e=.07;
t=50;
f=@(E) (w.*t)+(e.*sin(E))-E;
a=1.5.*10.^8;
b=-1.*10.^8;
error=10.^-3;
E_bi_a=bisectionMethod(f,a,b,error);
disp('The value for E for the bisection scheme is')
disp(E_bi_a)
end
%% 2.b. Implement the Newton-Raphson scheme to solve for E.
for i=1:1
w=.017; %set up values
e=.07;
t=50;
f=@(E) (w.*t)+(e.*sin(E))-E;
delta_t=.001;
error=10.^-3;
x_guess=1;
rootval=f(E_bi_a); %based on previous solution as found by bisection scheme
E_r=newton_ralph(f,delta_t,error,x_guess,rootval);
disp('The value for E for the bisection scheme is')
disp(E_r)
end
%% 2.c. For a range of times t = [0, 100] days and the parameters...
for i=1:1
for x=1:100 %Bisection Scheme
w=.017; %set up values
e=.07;
t=x;
f=@(Bisection_E) (w.*t)+(e.*sin(Bisection_E))-Bisection_E;
a=1.5.*10.^8;
b=-1.*10.^8;
error=10.^-3;
Bisection_E(x)=bisectionMethod(f,a,b,error);
end 
for x=1:100 %Newton-Ralphson Scheme
w=.017; %set up values
e=.07;
t=x;
f=@(NR_E) (w.*t)+(e.*sin(NR_E))-NR_E;
delta_t=.001;
error=10.^-3;
x_guess=1;
rootval=6.2381e-04; %based on previous solution as found by bisection scheme
NR_E(x)=newton_ralph(f,delta_t,error,x_guess,rootval);
end 
a=1.5.*10.^8;
e=.07;
w=.017;
for x=1:100
b_E=Bisection_E(x);
n_E=NR_E(x);
x_b(x)= a.*(cos(b_E)-e);
y_b(x)= a.*sqrt(1-(e.^2)).*sin(b_E);
x_n(x)= a.*(cos(n_E)-e);
y_n(x)= a.*sqrt(1-(e.^2)).*sin(n_E);
end
figure(1)
subplot(2,1,1)
plot(x_b,y_b,'b',-x_b,y_b,'b',-x_b,-y_b,'b',x_b,-y_b,'b')
title({'Bisection Scheme Derived Orbit'});
xlabel('x coordinates')
ylabel('y coordinates')
subplot(2,1,2)
plot(x_n,y_n,'r',-x_n,y_n,'r',-x_n,-y_n,'r',x_n,-y_n,'r')
title({'Newton-Ralphson Scheme Derived Orbit'});
xlabel('x coordinates')
ylabel('y coordinates')
figure(2)
plot(x_b,y_b,'b',x_n,y_n,'r')
title({'Bisection and Newton-Ralphson Scheme Derived Orbits, First Quadrant Only'});
xlabel('x coordinates')
ylabel('y coordinates')

for i=1:1 %figure 2 and argument for nondifferentiation betweem NR & B schemes
for x=1:100 %Bisection Scheme
w=.017; %set up values
e=.07;
t=x;
f=@(Bisection_E) (w.*t)+(e.*sin(Bisection_E))-Bisection_E;
a=1.5.*10.^8;
b=-1.*10.^8;
error=10.^-3;
Bisection_E(x)=bisectionMethod(f,a,b,error);
end 
for x=1:100 %Newton-Ralphson Scheme
w=.017; %set up values
e=.07;
t=x;
f=@(NR_E) (w.*t)+(e.*sin(NR_E))-NR_E;
delta_t=.001;
error=10.^-3;
x_guess=1;
rootval=3; %based on previous solution as found by bisection scheme
NR_E(x)=newton_ralph(f,delta_t,error,x_guess,rootval);
end 
a=1.5.*10.^8;
e=.07;
w=.017;
for x=1:100
b_E=Bisection_E(x);
n_E=NR_E(x);
x_b(x)= a.*(cos(b_E)-e);
y_b(x)= a.*sqrt(1-(e.^2)).*sin(b_E);
x_n(x)= a.*(cos(n_E)-e);
y_n(x)= a.*sqrt(1-(e.^2)).*sin(n_E);
end
figure(3)
plot(x_b,y_b,'b',-x_b,y_b,'b',-x_b,-y_b,'b',x_b,-y_b,'b')
hold on
plot(x_n,y_n,'r',-x_n,y_n,'r',-x_n,-y_n,'r',x_n,-y_n,'r')
title({'Newton-Ralphson Scheme Derived Orbit changes with different root val'});
xlabel('x coordinates')
ylabel('y coordinates')
legend({'Bisection-Scheme Derived Orbit','Newton-Ralphson Scheme Derived Orbit'})
    
disp('Because the rootval part of the newton-ralphson scheme was derived from')
disp('the answer for the bisection scheme, the values for E that I generated are')
disp('very close. If I change the value of rootval in the Newton-Ralphson')
disp('scheme, the results change dramatically, as seen in figure 3. I also')
disp('included the plots as not reflected over the other 3 quadrants, because')
disp('although this makes less sense in terms of the actual shape of an orbit,')
disp('this is what the functions originally outputted. You can''t see the ')
disp('Newton-Ralphson scheme when I do this- that''s because I used the original')
disp('root value that I derived from the Bisection Scheme for the NR scheme, like in figure 1).')
disp('I realized later that running x for longer creates the full orbit that''s')
disp('physically reasonable, but I thought I''d keep these figures to show my')
disp('thought process.')
end 
end
