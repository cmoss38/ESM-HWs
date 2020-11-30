%% 1. Nuclear Decay Equation
for i=1:1
%%a. Implement the Backward Euler method to numerically solve...
for i=1:1
%% n= 20
for i=1:1
n=20;
N_0=1;
t_s=0;
t_f=20;
lambda=[.1,1.1,2.1];

N=zeros(n+1,1);
N(1)=N_0;
N_exact=zeros(n+1,1);
N_exact(1)=N_0;

t=zeros(n+1,1);
error=zeros(n+1,1);
delta_t=(t_f-t_s)/n;

for j=1:length(lambda)
    l=lambda(j);
for i=1:n
    t(i+1)=t(i)+delta_t;
    N(i+1)=N(i)/(1+(l*delta_t));
    N_exact(i+1)=N_0*exp(-l*t(i+1));
    error(i+1)=error(i)+abs(N_exact(i+1)-N(i+1));
    
    figure(1)
    subplot(6,3,j)
    
    plot(t,N,'ko-','linewidth',3)
    hold on
    plot(t,N_exact,'r-','linewidth',3)
    xlabel('time')
    ylabel('N(t)')
    k=sprintf('Backward Euler Method for Decay Equation, delt=1, lambda= %.1f',l);
    title(k)
    legend({'Numerical Solution','Analytical Solution'})
    hold off
   
end

end
end
%% n=10
for i=1:1
n=10;
N_0=1;
t_s=0;
t_f=20;
lambda=[.1,1.1,2.1];

N=zeros(n+1,1);
N(1)=N_0;
N_exact=zeros(n+1,1);
N_exact(1)=N_0;

t=zeros(n+1,1);
error=zeros(n+1,1);
delta_t=(t_f-t_s)/n;
figs2=[4,5,6];

for j=1:length(lambda)
    l=lambda(j);
    ll=figs2(j);
for i=1:n
    t(i+1)=t(i)+delta_t;
    N(i+1)=N(i)/(1+(l*delta_t));
    N_exact(i+1)=N_0*exp(-l*t(i+1));
    error(i+1)=error(i)+abs(N_exact(i+1)-N(i+1));
    
    %figure(ll)
    subplot(6,3,ll)
    plot(t,N,'ko-','linewidth',3)
    hold on
    plot(t,N_exact,'g-','linewidth',3)
    xlabel('time')
    ylabel('N(t)')
    k=sprintf('Backward Euler Method for Decay Equation, delt=2, lambda= %.1f',l);
    title(k)
    legend({'Numerical Solution','Analytical Solution'})
    hold off
end
end
end
end
%%b. Implement the Centered Euler method to numerically solve... &
%%compare to the Backwards Euler method
for i=1:1
%% n= 20
for i=1:1
n=20;
N_0=1;
t_s=0;
t_f=20;
lambda=[.1,1.1,2.1];

N=zeros(n+2,1);
N(1)=N_0;
N_exact=zeros(n+2,1);
N_exact(1)=N_0;

t=zeros(n+1,1);
error=zeros(n+1,1);
delta_t=(t_f-t_s)/n;
figs=[7,8,9];

for j=1:length(lambda)
    l=lambda(j);
for i=1:n
    t(i+1)=t(i)+delta_t;
    N(i+2)=N(i+1)-N(i)/((2*delta_t));
    N_exact(i+1)=N_0*exp(-l*t(i+1));
    error(i+1)=error(i)+abs(N_exact(i+1)-N(i+1));
    
%     figure(figs(j))
    subplot(6,3,figs(j))
    plot(t,N(1:end-1),'ko-','linewidth',3)
    hold on
    plot(t,N_exact(1:end-1),'c-','linewidth',3)
    xlabel('time')
    ylabel('N(t)')
    k=sprintf('Centered Euler Method for Decay Equation, delt=1, lambda= %.1f',l);
    title(k)
    legend({'Numerical Solution','Analytical Solution'})
    hold off
end
end
end
%% n=10
for i=1:1
n=10;
N_0=1;
t_s=0;
t_f=20;
lambda=[.1,1.1,2.1];

N=zeros(n+2,1);
N(1)=N_0;
N_exact=zeros(n+2,1);
N_exact(1)=N_0;

t=zeros(n+1,1);
error=zeros(n+1,1);
delta_t=(t_f-t_s)/n;
figs=[10,11,12];
for j=1:length(lambda)
    l=lambda(j);
for i=1:n
    t(i+1)=t(i)+delta_t;
    N(i+2)=N(i+1)-N(i)/((2*delta_t));
    N_exact(i+1)=N_0*exp(-l*t(i+1));
    error(i+1)=error(i)+abs(N_exact(i+1)-N(i+1));
    
%     figure(figs(j))
    subplot(6,3,figs(j))
    plot(t,N(1:end-1),'ko-','linewidth',3)
    hold on
    plot(t,N_exact(1:end-1),'y-','linewidth',3)
    xlabel('time')
    ylabel('N(t)')
    k=sprintf('Centered Euler Method for Decay Equation, delt=2, lambda= %.1f',l);
    title(k)
    legend({'Numerical Solution','Analytical Solution'})
    hold off
end
end
end
disp('The Centered Euler Scheme is much less stable than the Backward Euler Scheme for all del(t) and lambda values.')
end
%%c. Calculate the stability condition for the Forward, Backward and Centered Euler methods
for i=1:1
%% n=20
for i=1:1
n=20;
N_0=1;
t_s=0;
t_f=20;
lambda=[.1,1.1,2.1];

N=zeros(n+1,1);
N(1)=N_0;
N_exact=zeros(n+1,1);
N_exact(1)=N_0;

t=zeros(n+1,1);
error=zeros(n+1,1);
delta_t=(t_f-t_s)/n;
figs=[13,14,15];

for j=1:length(lambda)
    l=lambda(j);
for i=1:n
    t(i+1)=t(i)+delta_t;
    N(i+1)=N(i)+(-l*(N(i))*delta_t);
    N_exact(i+1)=N_0*exp(-l*t(i+1));
    error(i+1)=error(i)+abs(N_exact(i+1)-N(i+1));
    
    
%     figure(figs(j))
    subplot(6,3,figs(j))
    plot(t,N,'ko-','linewidth',3)
    hold on
    plot(t,N_exact,'p-','linewidth',3)
    xlabel('time')
    ylabel('N(t)')
    k=sprintf('Forward Euler Method for Decay Equation, delt=1, lambda= %.1f',l);
    title(k)
    legend({'Numerical Solution','Analytical Solution'})
    hold off
end
end
end
%% n=10
for i=1:1
n=10;
N_0=1;
t_s=0;
t_f=20;
lambda=[.1,1.1,2.1];

N=zeros(n+1,1);
N(1)=N_0;
N_exact=zeros(n+1,1);
N_exact(1)=N_0;

t=zeros(n+1,1);
error=zeros(n+1,1);
delta_t=(t_f-t_s)/n;
figs=[16,17,18];

for j=1:length(lambda)
    l=lambda(j);
for i=1:n
    t(i+1)=t(i)+delta_t;
    N(i+1)=N(i)+(-l*(N(i))*delta_t);
    N_exact(i+1)=N_0*exp(-l*t(i+1));
    error(i+1)=error(i)+abs(N_exact(i+1)-N(i+1));
    
    %figure(figs(j))
    subplot(6,3,figs(j))
    plot(t,N,'ko-','linewidth',3)
    hold on
    plot(t,N_exact,'bl-','linewidth',3)
    xlabel('time')
    ylabel('N(t)')
    k=sprintf('Forward Euler Method for Decay Equation, delt=2, lambda= %.1f',l);
    title(k)
    legend({'Numerical Solution','Analytical Solution'})
    hold off
    
end
end
end
end
disp(' ')
disp('All schemes are unconditionally stable based on the parameters inputted,')
disp('EXCEPT for several of the Forward Euler Method scenarios, making the')
disp('Forward Euler method conditionally stable for the Nuclear Decay')
disp('Equation.')
disp('Forward Euler method for del(t)=1, lambda=2.1; del(t)=2, lambda=1.1;')
disp('del(t)=2, lambda=2.1 are unstable.')
disp(' ')
disp('Numerical instability causes the graphs to become more jagged an approach')
disp('infinity rather than zero.')
disp(' ')
disp('I would use Forward Euler Scheme to graph this because the mumerical')
disp('solution follows the analytical solution the closest across all parameter')
disp('adjustments.')
end

%% 2. Lorenz Equations/Attractor
for i=1:1
%%a. Implement the Predictor-Corrector method to numerically solve
for i=1:1
t(1)=0;  
x(1)=1;
y(1)=1;
z(1)=1;
sigma=10;  
rho=28;
beta=8/3;
delta_t=0.01;   
t=0:delta_t:100;
f=@(t,x,y,z) sigma*(y-x);   
g=@(t,x,y,z) x*rho-x.*z-y;
p=@(t,x,y,z) x.*y-beta*z;
for i=1:(length(t)-1)
    
    k1=f(t(i),x(i),y(i),z(i));
    l1=g(t(i),x(i),y(i),z(i));
    m1=p(t(i),x(i),y(i),z(i));
    
    k2=f(t(i)+delta_t/2,((x(i)+k1)*delta_t),((y(i)+l1)*delta_t),((z(i)+m1)*delta_t));     
    l2=g(t(i)+delta_t/2,((x(i)+k1)*delta_t),((y(i)+l1)*delta_t),((z(i)+m1)*delta_t));
    m2=p(t(i)+delta_t/2,((x(i)+k1)*delta_t),((y(i)+l1)*delta_t),((z(i)+m1)*delta_t));
      
    x(i+1)=x(i)+delta_t*(k1+k2)/2;
    y(i+1)=y(i)+delta_t*(l1+l2)/2;
    z(i+1)=z(i)+delta_t*(m1+m2)/2;
end
figure(2)
plot3(x,y,z,'b')
hold on

end
%b. Implement the Predictor-Corrector method to numerically solve, (x,y,z)=(1.000000001,1,1)

for i=1:1
t(1)=0;  
x(1)=1.000000001;
y(1)=1;
z(1)=1;
sigma=10;  
rho=28;
beta=8/3;
delta_t=0.01;   
t=0:delta_t:100;
f=@(t,x,y,z) sigma*(y-x);   
g=@(t,x,y,z) x*rho-x.*z-y;
p=@(t,x,y,z) x.*y-beta*z;
for i=1:(length(t)-1)
    
    k1=f(t(i),x(i),y(i),z(i));
    l1=g(t(i),x(i),y(i),z(i));
    m1=p(t(i),x(i),y(i),z(i));
    
    k2=f(t(i)+delta_t/2,((x(i)+k1)*delta_t),((y(i)+l1)*delta_t),((z(i)+m1)*delta_t));     
    l2=g(t(i)+delta_t/2,((x(i)+k1)*delta_t),((y(i)+l1)*delta_t),((z(i)+m1)*delta_t));
    m2=p(t(i)+delta_t/2,((x(i)+k1)*delta_t),((y(i)+l1)*delta_t),((z(i)+m1)*delta_t));
      
    x(i+1)=x(i)+delta_t*(k1+k2)/2;
    y(i+1)=y(i)+delta_t*(l1+l2)/2;
    z(i+1)=z(i)+delta_t*(m1+m2)/2;
end

plot3(x,y,z,'r')

end
%%c. Implement the Fourth-Order Runge-Kutta (RK4) method to numerically solve
for i=1:1
t(1)=0;  
x(1)=1;
y(1)=1;
z(1)=1;
sigma=10;  
rho=28;
beta=8/3;
delta_t=0.01;   
t=0:delta_t:100;
f=@(t,x,y,z) sigma*(y-x);   
g=@(t,x,y,z) x*rho-x.*z-y;
p=@(t,x,y,z) x.*y-beta*z;
for i=1:(length(t)-1)
    
    k1=f(t(i),x(i),y(i),z(i));
    l1=g(t(i),x(i),y(i),z(i));
    m1=p(t(i),x(i),y(i),z(i));
    
    k2=f(t(i)+delta_t/2,(x(i)+0.5*k1*delta_t),(y(i)+(0.5*l1*delta_t)),(z(i)+(0.5*m1*delta_t)));     
    l2=g(t(i)+delta_t/2,(x(i)+0.5*k1*delta_t),(y(i)+(0.5*l1*delta_t)),(z(i)+(0.5*m1*delta_t)));
    m2=p(t(i)+delta_t/2,(x(i)+0.5*k1*delta_t),(y(i)+(0.5*l1*delta_t)),(z(i)+(0.5*m1*delta_t)));
    k3=f(t(i)+delta_t/2,(x(i)+0.5*k2*delta_t),(y(i)+(0.5*l2*delta_t)),(z(i)+(0.5*m2*delta_t)));
    l3=g(t(i)+delta_t/2,(x(i)+0.5*k2*delta_t),(y(i)+(0.5*l2*delta_t)),(z(i)+(0.5*m2*delta_t)));
    m3=p(t(i)+delta_t/2,(x(i)+0.5*k2*delta_t),(y(i)+(0.5*l2*delta_t)),(z(i)+(0.5*m2*delta_t)));
    k4=f(t(i)+delta_t,(x(i)+k3*delta_t),(y(i)+l3*delta_t),(z(i)+m3*delta_t));
    l4=g(t(i)+delta_t,(x(i)+k3*delta_t),(y(i)+l3*delta_t),(z(i)+m3*delta_t));
    m4=p(t(i)+delta_t,(x(i)+k3*delta_t),(y(i)+l3*delta_t),(z(i)+m3*delta_t));
      
    x(i+1)=x(i)+delta_t*(k1+2*k2+2*k3+k4)/6;
    y(i+1)=y(i)+delta_t*(l1+2*l2+2*l3+l4)/6;
    z(i+1)=z(i)+delta_t*(m1+2*m2+2*m3+m4)/6;
end

plot3(x,y,z,'g')
title('Lorenz Attractor Solved Numerically')
xlabel('x')
ylabel('y')
zlabel('z')
legend({'Predictor-Corrector Method, (x1,y1,z1)=(1,1,1)','Predictor-Corrector Method, (x1,y1,z1)=(1.000000001,1,1)','Runge-Kutta Method, (x1,y1,z1)=(1,1,1)'})
hold off 
disp(' ')
disp('Even the *slightest* perturbation or change in method will change these plots')
end
end