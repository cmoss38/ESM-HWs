for i = 1:1 %1.Consider the Integral... 
%% a. Plot the function f(x)
x=linspace(1,100); %just to get some values to run the fxn through
f=2.*exp(x)+x.^4-cos(x); %fxn given in problem

figure(1)
plot(x,f)
hold on
xlabel('x values')
ylabel('y values')
title('f(x) as integrated by MATLAB')
%% b. Find the exact value of the integral
x=0; %value to plug into the integral of f
F=(1./5).*x.^5+2.*exp(x)-sin(x); %integrated f to get F
disp('The exact integral is (1/5)(pi^5)+2e^pi-sin(pi)-2');

%% c. Write a Matlab script that implements the Midpoint Rule for numerical integration of the function f(x)
% I. set up 
b=0; %lower bound, adjustable
a=pi; %upper bound, adjustable
n=1:50; %number of subintervals for each iteration of the Midpoint rule
I=[]; %initializing the loop that we're going to use to plot against the exact value of I

% II. loops for each plot
for j=1:length(n) %this loop runs through all the values we've put into the n variable
    deltax=(a-b)/n(j); %creating integral using midpoint rule 
    delta=b:deltax:a; %^^^
    %F= @(x) (1./5).*x.^5+2.*exp(x)-sin(x); % learned how to make a fxn handle https://www.mathworks.com/matlabcentral/answers/91250-midpoint-rule-for-integration
    f= @(x) 2.*exp(x)+x.^4-cos(x);
    mpr=0; %initializing midpoint rule  
    for k=1:n(j)
        x_k=b+(2.*(k-1)+1)./2.*deltax;
        mpr=mpr+(f(x_k).*deltax); %adding up values for each iteration of mpr
    end
    I=[I mpr]; %compiles every iteration together
end

% III. comparing differences 
I_exact=(1./5).*(pi.^5)+2.*exp(pi)-sin(pi)-2; %value we found earlier but we're actually going to let matlab run it this time
diff=abs(I-I_exact);

figure(2)
plot(n,diff) 
xlabel('Number of Time Steps')
ylabel('Total Error')
title('Error from Midpoint Rule Decreases with More Time Steps')
mpr_total=sum(diff);
%% d. Write a Matlab script that implements the Trapezoid Rule for numerical integration of the function f(x)
% I. set up 
b=0; %lower bound, adjustable
a=pi; %upper bound, adjustable
n=1:50; %subintervals for each iteration of the Midpoint rule
I=[]; %initializing the loop that we're going to use to plot against the exact value of I
% II. loops for each plot
for j=1:length(n) %this loop runs through all the values we've put into the n variable
deltax=(a-b)/n(j); %creating integral using midpoint rule 
x=b:deltax:a; %^^^
f= @(x) 2.*exp(x)+x.^4-cos(x);
tr=0; %initializing midpoint rule  
    for k=1:n(j)
        x_k=b+(k-1).*deltax;
        x_k1=b+(k).*deltax;
        tr=tr+(f(x_k1)+f(x_k))/2.*deltax; %adding up values for each iteration of mpr
    end
I=[I tr]; %compiles every iteration together
end

% III. comparing differences 
I_exact=(1./5).*(pi.^5)+2.*exp(pi)-sin(pi)-2; %value we found earlier but we're actually going to let matlab run it this time
diff=abs(I-I_exact);

figure(3)
plot(n,diff)
xlabel('Number of Time Steps')
ylabel('Total Error')
title('Error from Trapezoid Rule Decreases with More Time Steps')
tr_total=sum(diff);

%% e. Compare how well these two methods numerically approximate the exact value of this integral
if tr_total < mpr_total
    disp('The trapezoid rule approximates the value for the integral better which is what I expected.')
end
if tr_total > mpr_total
    disp('The midpoint rule approximates the value for the integral better which is what I expected.')
end
end

clear all
close all
 
for i=1:1 %2. Given an instantaneous change of mass distribution...
%% a. Calculate the expected sea level change
p_e=5.5; %avg density of earth, g/cm^3
p_w=1; %density of ocean water, g/cm^3
b=pi/2; %upper bound, adjustable
a=-pi/2; %lower bound, adjustable
n=500; %subintervals for each iteration of the Midpoint rule
I=[]; %initializing the loop that we're going to use to plot against the exact value of I
deltax=(b-a)/n; %creating integral using midpoint rule 
x=a:deltax:b;
theta=.59; %as given in problem

pt1=(-1./(2.*sin(abs(x-theta)./2)+10.^-3))-6.5;
pt2=exp((-10.*(x+(pi./2)).^2)).*deltax;
f= @(x) pt1.*pt2;

tr=0; %initializing trapezoid rule  
    for k=1:n
        x_k=a+(k-1).*deltax;
        x_k1=a+(k).*deltax;
        tr=tr+(f(x_k1)+f(x_k))/2.*deltax; %adding up values for each iteration of tr
    end
I=[I tr]; %compiles every iteration together
S=sum(I(2:end)); %get rid of extra value at beginnning
disp('the value of S is -6.18') 
%% b. calculate and plot the expected sea level change (S) as a function of latitude over 100 evenly spaced values of 0, in the range [-1:5;1:5].
p_e=5.5; %avg density of earth, g/cm^3
p_w=1; %density of ocean water, g/cm^3
b=pi/2; %upper bound, adjustable
a=-pi/2; %lower bound, adjustable
n=500; %subintervals for each iteration of the Midpoint rule
I=[]; %initializing the loop that we're going to use to plot against the exact value of I
deltax=(b-a)/n; %creating integral using midpoint rule 
x=a:deltax:b;
theta=linspace(-1.5,1.5); %as given in problem

for x_t=1:length(theta)
pt1=(-1./(2.*sin(abs(x-theta(x_t))./2)+10.^-3))-6.5;
pt2=exp((-10.*(x+(pi./2)).^2)).*deltax;
f= @(x) pt1.*pt2;

tr=0; %initializing trapezoid rule  
    for k=1:n
        x_k=a+(k-1).*deltax;
        x_k1=a+(k).*deltax;
        tr=tr+(f(x_k1)+f(x_k))/2.*deltax; %adding up values for each iteration of tr
    end
I=[I;tr]; %compiles every iteration together
d=sum(I');
r=sum(I);
conv_lat=.59./33.7490;
end

figure(4) %right answer
plot(theta/conv_lat,d)
xlabel('Latitude')
ylabel('Sea Level Change (cm)')
title('Sea level as a Function of Latitude')
%% c. Explain why the sea level change curve in part (b) has the sign and shape that it does
disp('The sea level change curve has the sign and shape that it does')
disp('because the addition of mass (i.e. ice from snow accumulation) causes isostatic depression at the South Pole, which causes the sea around it to')
disp('sink in elevation. The sea farther away from it isn''t close enough for the elevation effect of isostatic rebound at a regional level, so sea level change is less')
disp('negative, although it''s still negative (because water from the ocean directly added to ice mass on Antarctica''s ice sheet.')
end