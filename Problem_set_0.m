%% ESM HW 0             Hudson Moss 

clear all
close all
%% 1
for x=1
%1.a. nonvectorized version
k=300; %value given in problem 
s=0;
n=365;
for d=1:n %vector of every day in a year
    s_end=(cos((2*pi*d)/(365))+2)/k; %equation given in problem
    s=s+s_end; %compiles the snowfall for each day
end
disp(s);

%1.b. vectorized version
k=300; 
s=0;
d=1:365;
s2=sum(cos((2.*pi.*d)./(365))+2)./k; %same as previous but less computationally expensive
disp(s2);

%1.c. climate controlled 
d=1:365; %number of days in a year
s=3; %initializing the accumulation
k=300; %initializing the climate factor 
while s>2  %we are trying to find the maximum value of k for which s is less than 2 m.
        s=sum((cos((2.*pi.*d)./(365))+2)./k);
        k=k+1;
end
disp(k)

clear all
close all
end
%% 2
for x=1

%2.a.fibonacci values 20-200 
F=fibonacci(0:20);
F=F(F>=20 & F<=200); %figures out which values don't fit in range and deletes them

%2.b.fibo.dat
a_n=fibonacci(0:19);
n=1:20;
fib_file=[n ;a_n]; %gives index number and actual value starting at 0
fileName = fopen('fibo.dat','w'); %this writes to a file
fprintf(fileName, '%f %f\n', fib_file);
fclose(fileName)
type fibo.dat %spits out values in file for me to check them

clear all
close all 
end
%% 3 
for x=1
%3.a. matrix operations
for n=1
vect=[1 1 1 4;1 -2 3 -6;2 3 1 7]; %this is the original set of values given in the problem
%1
vect(2,:)=vect(2,:)+(vect(1,:).*-1); %adding -1 times the first row to the second row
%2
vect(3,:)=vect(3,:)+(vect(1,:).*-2); %adding -2 times the first row to the third row
%3
vect(2,:)=vect(2,:).*(-1/3); %dividing the second row by -3
%4
vect(3,:)=vect(3,:)+(vect(2,:).*-1); %adding -1 times the second row to the third 
%5
vect(3,:)=vect(3,:).*(3/13); %multipling the third row by 3/13
%6
vect(2,:)=vect(2,:)+(vect(3,:).*(4/3)); %adding 4/3 times the third row to the second
%7
vect(1,:)=vect(1,:)+vect(3,:); %adding the third row to the first 
%8
vect(1,:)=vect(1,:)+(vect(2,:).*-1); %adding -1 times the second row to the first 

x=vect(1,4); %answers found from reducing matrix to reduced row echelon form
y=vect(2,4);
z=vect(3,4);
end %I manually did the operations to put it in RREF but there's a better way below
%or literally do it matrix operations? instructions unclear whoops im dumb
a=[1 1 -1;1 -2 3; 2 3 1]; %x,y,z coefficients for each equation, respectively given
b=[4 -6 7]'; %solution sets for each equation
n=inv(a)*b; %finds solution 
x=n(1,:) %these just spit out the respective x,y,z values for the system
y=n(2,:)
z=n(3,:)

clear all
close all
end
%% 4 
for x=1
%4.a.
y=[0:.1:3];
z=.25.*y+.5.*y.*exp(-.5^2-y.^2); %equation for surface height 
figure(1)
plot(y,z)
title('Surface Height as a Function of Y direction')
ylabel('Surface Height (m)')
xlabel('Y direction (m)')

%4.b.
figure(2)
e=2.72;
[x,y] = meshgrid(0:0.1:2,0:0.1:2); %this just makes a bunch of points that I'll run z through
z=.25.*(y)+(x).*(y).*(e.^((-x.^2)-y.^2));
surf(x,y,z,'edgecolor','none')
hold on 
alpha .8 %this is just to make the plot transparent

%4.c. 
[val,idx]=max(z(:)) ;  % get maximum in Z 
x=x(idx) ;  % get x for zmax 
y=y(idx) ;  % get y for zmax 
hold on
plot3(x,y,val,'o','markersize',20) ;
title('Surface Height as a Function of X and Y, Max Z Value Circled')
ylabel('X direction (m)')
xlabel('Y direction (m)')
zlabel('Surface Height (m)')

end
    