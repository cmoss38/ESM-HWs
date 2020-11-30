%% HW 5 Finally Done!!
% 2. Write a MATLAB script that solves the advection equation by applying the Forward Euler upwind scheme
%% c=.1
for i=1:1
    %t=100000
%I. Set up
u=0.001;    
c=0.1;      %courant number

delta_x=2;
x=0:delta_x:2000;
delta_t=(c*delta_x)/u;
t=0:delta_t:2589408; %number of sec/30 days
J_x=zeros(1000,1);

%II. Matrix
for i=1:length(x)-1
    J_x(i)= exp(-0.0001.*(x(i)-300)^2);
    i=i+1;
end

A=ones(1000,1);
JM=spdiags([c.*A (1-c)*A],-1:0,1000,1000);
full(JM);

t_final=100000;
tnew=0:delta_t:t_final;

for k= 1:length(tnew)
    Jnew=JM*J_x;
    J_x=Jnew;
    k=k+1;
end

%III. Plot
x=0:delta_x:1999;
figure(1)
plot(x,J_x)
hold on

%==========================================================================
%t=200000
%I. Set up
u=0.001;    
c=0.1;      %courant number

delta_x=2;
x=0:delta_x:2000;
delta_t=(c*delta_x)/u;
t=0:delta_t:2589408; %number of sec/30 days
J_x=zeros(1000,1);

%II. Matrix
for i=1:length(x)-1
    J_x(i)= exp(-0.0001.*(x(i)-300)^2);
    i=i+1;
end

A=ones(1000,1);
JM=spdiags([c.*A (1-c)*A],-1:0,1000,1000);
full(JM);

t_final=200000;
tnew=0:delta_t:t_final;

for k= 1:length(tnew)
    Jnew=JM*J_x;
    J_x=Jnew;
    k=k+1;
end

%III. Plot
x=0:delta_x:1999;
plot(x,J_x)
hold on

%==========================================================================
%t=300000
%I. Set up
u=0.001;    
c=0.1;      %courant number

delta_x=2;
x=0:delta_x:2000;
delta_t=(c*delta_x)/u;
t=0:delta_t:2589408; %number of sec/30 days
J_x=zeros(1000,1);

%II. Matrix
for i=1:length(x)-1
    J_x(i)= exp(-0.0001.*(x(i)-300)^2);
    i=i+1;
end

A=ones(1000,1);
JM=spdiags([c.*A (1-c)*A],-1:0,1000,1000);
full(JM);

t_final=300000;
tnew=0:delta_t:t_final;

for k= 1:length(tnew)
    Jnew=JM*J_x;
    J_x=Jnew;
    k=k+1;
end

%III. Plot
x=0:delta_x:1999;
plot(x,J_x)
hold on
%==========================================================================
%t=400000
%I. Set up
u=0.001;    
c=0.1;      %courant number

delta_x=2;
x=0:delta_x:2000;
delta_t=(c*delta_x)/u;
t=0:delta_t:2589408; %number of sec/30 days
J_x=zeros(1000,1);

%II. Matrix
for i=1:length(x)-1
    J_x(i)= exp(-0.0001.*(x(i)-300)^2);
    i=i+1;
end

A=ones(1000,1);
JM=spdiags([c.*A (1-c)*A],-1:0,1000,1000);
full(JM);

t_final=400000;
tnew=0:delta_t:t_final;

for k= 1:length(tnew)
    Jnew=JM*J_x;
    J_x=Jnew;
    k=k+1;
end

%III. Plot
x=0:delta_x:1999;
plot(x,J_x)
hold on
%==========================================================================
%t=500000
%I. Set up
u=0.001;    
c=0.1;      %courant number

delta_x=2;
x=0:delta_x:2000;
delta_t=(c*delta_x)/u;
t=0:delta_t:2589408; %number of sec/30 days
J_x=zeros(1000,1);

%II. Matrix
for i=1:length(x)-1
    J_x(i)= exp(-0.0001.*(x(i)-300)^2);
    i=i+1;
end

A=ones(1000,1);
JM=spdiags([c.*A (1-c)*A],-1:0,1000,1000);
full(JM);

t_final=500000;
tnew=0:delta_t:t_final;

for k= 1:length(tnew)
    Jnew=JM*J_x;
    J_x=Jnew;
    k=k+1;
end

%III. Plot
x=0:delta_x:1999;
plot(x,J_x)
hold on
xlabel('Distance (km)')
ylabel('Concentration (normalized)')
title('Concentration over time, C=.1')
legend({'t=100000','t=200000','t=300000','t=400000','t=500000'})
    
end
%% c=.5
for i=1:1
    %t=100000
%I. Set up
u=0.001;    
c=0.5;      %courant number

delta_x=2;
x=0:delta_x:2000;
delta_t=(c*delta_x)/u;
t=0:delta_t:2589408; %number of sec/30 days
J_x=zeros(1000,1);

%II. Matrix
for i=1:length(x)-1
    J_x(i)= exp(-0.0001.*(x(i)-300)^2);
    i=i+1;
end

A=ones(1000,1);
JM=spdiags([c.*A (1-c)*A],-1:0,1000,1000);
full(JM);

t_final=100000;
tnew=0:delta_t:t_final;

for k= 1:length(tnew)
    Jnew=JM*J_x;
    J_x=Jnew;
    k=k+1;
end

%III. Plot
x=0:delta_x:1999;
figure(2)
plot(x,J_x)
hold on

%==========================================================================
%t=200000
%I. Set up
u=0.001;    

delta_x=2;
x=0:delta_x:2000;
delta_t=(c*delta_x)/u;
t=0:delta_t:2589408; %number of sec/30 days
J_x=zeros(1000,1);

%II. Matrix
for i=1:length(x)-1
    J_x(i)= exp(-0.0001.*(x(i)-300)^2);
    i=i+1;
end

A=ones(1000,1);
JM=spdiags([c.*A (1-c)*A],-1:0,1000,1000);
full(JM);

t_final=200000;
tnew=0:delta_t:t_final;

for k= 1:length(tnew)
    Jnew=JM*J_x;
    J_x=Jnew;
    k=k+1;
end

%III. Plot
x=0:delta_x:1999;
plot(x,J_x)
hold on

%==========================================================================
%t=300000
%I. Set up
u=0.001;    

delta_x=2;
x=0:delta_x:2000;
delta_t=(c*delta_x)/u;
t=0:delta_t:2589408; %number of sec/30 days
J_x=zeros(1000,1);

%II. Matrix
for i=1:length(x)-1
    J_x(i)= exp(-0.0001.*(x(i)-300)^2);
    i=i+1;
end

A=ones(1000,1);
JM=spdiags([c.*A (1-c)*A],-1:0,1000,1000);
full(JM);

t_final=300000;
tnew=0:delta_t:t_final;

for k= 1:length(tnew)
    Jnew=JM*J_x;
    J_x=Jnew;
    k=k+1;
end

%III. Plot
x=0:delta_x:1999;
plot(x,J_x)
hold on
%==========================================================================
%t=400000
%I. Set up
u=0.001;    

delta_x=2;
x=0:delta_x:2000;
delta_t=(c*delta_x)/u;
t=0:delta_t:2589408; %number of sec/30 days
J_x=zeros(1000,1);

%II. Matrix
for i=1:length(x)-1
    J_x(i)= exp(-0.0001.*(x(i)-300)^2);
    i=i+1;
end

A=ones(1000,1);
JM=spdiags([c.*A (1-c)*A],-1:0,1000,1000);
full(JM);

t_final=400000;
tnew=0:delta_t:t_final;

for k= 1:length(tnew)
    Jnew=JM*J_x;
    J_x=Jnew;
    k=k+1;
end

%III. Plot
x=0:delta_x:1999;
plot(x,J_x)
hold on
%==========================================================================
%t=500000
%I. Set up
u=0.001;    

delta_x=2;
x=0:delta_x:2000;
delta_t=(c*delta_x)/u;
t=0:delta_t:2589408; %number of sec/30 days
J_x=zeros(1000,1);

%II. Matrix
for i=1:length(x)-1
    J_x(i)= exp(-0.0001.*(x(i)-300)^2);
    i=i+1;
end

A=ones(1000,1);
JM=spdiags([c.*A (1-c)*A],-1:0,1000,1000);
full(JM);

t_final=500000;
tnew=0:delta_t:t_final;

for k= 1:length(tnew)
    Jnew=JM*J_x;
    J_x=Jnew;
    k=k+1;
end

%III. Plot
x=0:delta_x:1999;
plot(x,J_x)
hold on
xlabel('Distance (km)')
ylabel('Concentration (normalized)')
title('Concentration over time, C=.5')
legend({'t=100000','t=200000','t=300000','t=400000','t=500000'})

end
%% c=1
for i=1:1
    %t=100000
%I. Set up
u=0.001;    
c=1;      %courant number

delta_x=2;
x=0:delta_x:2000;
delta_t=(c*delta_x)/u;
t=0:delta_t:2589408; %number of sec/30 days
J_x=zeros(1000,1);

%II. Matrix
for i=1:length(x)-1
    J_x(i)= exp(-0.0001.*(x(i)-300)^2);
    i=i+1;
end

A=ones(1000,1);
JM=spdiags([c.*A (1-c)*A],-1:0,1000,1000);
full(JM);

t_final=100000;
tnew=0:delta_t:t_final;

for k= 1:length(tnew)
    Jnew=JM*J_x;
    J_x=Jnew;
    k=k+1;
end

%III. Plot
x=0:delta_x:1999;
figure(3)
plot(x,J_x)
hold on

%==========================================================================
%t=200000
%I. Set up
u=0.001;    

delta_x=2;
x=0:delta_x:2000;
delta_t=(c*delta_x)/u;
t=0:delta_t:2589408; %number of sec/30 days
J_x=zeros(1000,1);

%II. Matrix
for i=1:length(x)-1
    J_x(i)= exp(-0.0001.*(x(i)-300)^2);
    i=i+1;
end

A=ones(1000,1);
JM=spdiags([c.*A (1-c)*A],-1:0,1000,1000);
full(JM);

t_final=200000;
tnew=0:delta_t:t_final;

for k= 1:length(tnew)
    Jnew=JM*J_x;
    J_x=Jnew;
    k=k+1;
end

%III. Plot
x=0:delta_x:1999;
plot(x,J_x)
hold on

%==========================================================================
%t=300000
%I. Set up
u=0.001;    

delta_x=2;
x=0:delta_x:2000;
delta_t=(c*delta_x)/u;
t=0:delta_t:2589408; %number of sec/30 days
J_x=zeros(1000,1);

%II. Matrix
for i=1:length(x)-1
    J_x(i)= exp(-0.0001.*(x(i)-300)^2);
    i=i+1;
end

A=ones(1000,1);
JM=spdiags([c.*A (1-c)*A],-1:0,1000,1000);
full(JM);

t_final=300000;
tnew=0:delta_t:t_final;

for k= 1:length(tnew)
    Jnew=JM*J_x;
    J_x=Jnew;
    k=k+1;
end

%III. Plot
x=0:delta_x:1999;
plot(x,J_x)
hold on
%==========================================================================
%t=400000
%I. Set up
u=0.001;    

delta_x=2;
x=0:delta_x:2000;
delta_t=(c*delta_x)/u;
t=0:delta_t:2589408; %number of sec/30 days
J_x=zeros(1000,1);

%II. Matrix
for i=1:length(x)-1
    J_x(i)= exp(-0.0001.*(x(i)-300)^2);
    i=i+1;
end

A=ones(1000,1);
JM=spdiags([c.*A (1-c)*A],-1:0,1000,1000);
full(JM);

t_final=400000;
tnew=0:delta_t:t_final;

for k= 1:length(tnew)
    Jnew=JM*J_x;
    J_x=Jnew;
    k=k+1;
end

%III. Plot
x=0:delta_x:1999;
plot(x,J_x)
hold on
%==========================================================================
%t=500000
%I. Set up
u=0.001;    

delta_x=2;
x=0:delta_x:2000;
delta_t=(c*delta_x)/u;
t=0:delta_t:2589408; %number of sec/30 days
J_x=zeros(1000,1);

%II. Matrix
for i=1:length(x)-1
    J_x(i)= exp(-0.0001.*(x(i)-300)^2);
    i=i+1;
end

A=ones(1000,1);
JM=spdiags([c.*A (1-c)*A],-1:0,1000,1000);
full(JM);

t_final=500000;
tnew=0:delta_t:t_final;

for k= 1:length(tnew)
    Jnew=JM*J_x;
    J_x=Jnew;
    k=k+1;
end

%III. Plot
x=0:delta_x:1999;
plot(x,J_x)
hold on
xlabel('Distance (km)')
ylabel('Concentration (normalized)')
title('Concentration over time, C=1')
legend({'t=100000','t=200000','t=300000','t=400000','t=500000'})

end
%% c=2
for i=1:1
    %t=100000
%I. Set up
u=0.001;    
c=2;      %courant number

delta_x=2;
x=0:delta_x:2000;
delta_t=(c*delta_x)/u;
t=0:delta_t:2589408; %number of sec/30 days
J_x=zeros(1000,1);

%II. Matrix
for i=1:length(x)-1
    J_x(i)= exp(-0.0001.*(x(i)-300)^2);
    i=i+1;
end

A=ones(1000,1);
JM=spdiags([c.*A (1-c)*A],-1:0,1000,1000);
full(JM);

t_final=100000;
tnew=0:delta_t:t_final;

for k= 1:length(tnew)
    Jnew=JM*J_x;
    J_x=Jnew;
    k=k+1;
end

%III. Plot
x=0:delta_x:1999;
figure(4)
plot(x,J_x,'linewidth',3)
hold on

%==========================================================================
%t=200000
%I. Set up
u=0.001;    

delta_x=2;
x=0:delta_x:2000;
delta_t=(c*delta_x)/u;
t=0:delta_t:2589408; %number of sec/30 days
J_x=zeros(1000,1);

%II. Matrix
for i=1:length(x)-1
    J_x(i)= exp(-0.0001.*(x(i)-300)^2);
    i=i+1;
end

A=ones(1000,1);
JM=spdiags([c.*A (1-c)*A],-1:0,1000,1000);
full(JM);

t_final=200000;
tnew=0:delta_t:t_final;

for k= 1:length(tnew)
    Jnew=JM*J_x;
    J_x=Jnew;
    k=k+1;
end

%III. Plot
x=0:delta_x:1999;
plot(x,J_x,'linewidth',3)
hold on

%==========================================================================
%t=300000
%I. Set up
u=0.001;    

delta_x=2;
x=0:delta_x:2000;
delta_t=(c*delta_x)/u;
t=0:delta_t:2589408; %number of sec/30 days
J_x=zeros(1000,1);

%II. Matrix
for i=1:length(x)-1
    J_x(i)= exp(-0.0001.*(x(i)-300)^2);
    i=i+1;
end

A=ones(1000,1);
JM=spdiags([c.*A (1-c)*A],-1:0,1000,1000);
full(JM);

t_final=300000;
tnew=0:delta_t:t_final;

for k= 1:length(tnew)
    Jnew=JM*J_x;
    J_x=Jnew;
    k=k+1;
end

%III. Plot
x=0:delta_x:1999;
plot(x,J_x,'linewidth',3)
hold on
%==========================================================================
%t=400000
%I. Set up
u=0.001;    

delta_x=2;
x=0:delta_x:2000;
delta_t=(c*delta_x)/u;
t=0:delta_t:2589408; %number of sec/30 days
J_x=zeros(1000,1);

%II. Matrix
for i=1:length(x)-1
    J_x(i)= exp(-0.0001.*(x(i)-300)^2);
    i=i+1;
end

A=ones(1000,1);
JM=spdiags([c.*A (1-c)*A],-1:0,1000,1000);
full(JM);

t_final=400000;
tnew=0:delta_t:t_final;

for k= 1:length(tnew)
    Jnew=JM*J_x;
    J_x=Jnew;
    k=k+1;
end

%III. Plot
x=0:delta_x:1999;
plot(x,J_x,'linewidth',3)
hold on
%==========================================================================
%t=500000
%I. Set up
u=0.001;    

delta_x=2;
x=0:delta_x:2000;
delta_t=(c*delta_x)/u;
t=0:delta_t:2589408; %number of sec/30 days
J_x=zeros(1000,1);

%II. Matrix
for i=1:length(x)-1
    J_x(i)= exp(-0.0001.*(x(i)-300)^2);
    i=i+1;
end

A=ones(1000,1);
JM=spdiags([c.*A (1-c)*A],-1:0,1000,1000);
full(JM);

t_final=500000;
tnew=0:delta_t:t_final;

for k= 1:length(tnew)
    Jnew=JM*J_x;
    J_x=Jnew;
    k=k+1;
end

%III. Plot
x=0:delta_x:1999;
plot(x,J_x)
hold on
xlabel('Distance (km)')
ylabel('Concentration (normalized)')
title('Concentration over time, C=2')
legend({'t=100000','t=200000','t=300000','t=400000','t=500000'})

end

disp('As courant number increases, you can see that Matlab is worse at refining')
disp('the spread of concentration of chemicals. As t increases for the c=.1 iteration of the code,')
disp('the curve for concentration over distance has increased in spread by the time')
disp('it hits t=500000. For the c=1 iteration, the concentration at the center of the')
disp('distribution has barely changed (because the resolution is much lower because of the')
disp('high courant number). For c=2, the graph is no longer physically coherent- c usually')
disp('has to be between 0 and 1 in order for the wave to make sense.')