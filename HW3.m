%% 1.b.
for i=1:1
%I. Setup
delta_t=.1; %yrs 
T_N(1)=260; %Kelvin
T_E(1)=260; %Kelvin
T_S(1)=260; %Kelvin
a_N=.7; %Northern Albedo
a_E=.5; %Equatorial Albedo
a_S=.7; %Southern Albedo
e=10E-7; %epsilon value 
S=1250; %flux from sun
D=.5; %factor for 
t=0:delta_t:1000; %generating time steps

for i=1:(length(t)-1)

    k1_T_N=((1-a_N).*S-e.*T_N(i).^4+D.*(T_E(i)-T_N(i))).*delta_t;
    k1_T_S=((1-a_S).*S-e.*T_S(i).^4+D.*(T_E(i)-T_S(i))).*delta_t;
    k1_T_E=((1-a_E).*S-e.*T_E(i).^4+D.*((2.*T_E(i))-T_N(i)-T_S(i))).*delta_t;
    
    k2_T_N=((1-a_N).*S-e.*(T_N(i)+k1_T_N)^4+D.*((T_E(i)+k1_T_E)-(T_N(i)-k1_T_N))).*delta_t;
    k2_T_S=((1-a_S).*S-e.*(T_S(i)+k1_T_S)^4+D.*((T_E(i)+k1_T_E)-(T_S(i)-k1_T_S))).*delta_t;
    k2_T_E=((1-a_E).*S-e.*(T_E(i)+k1_T_E)^4+D.*((2.*(T_E(i)+k1_T_E))-(T_N(i)-k1_T_N)-(T_S(i)-k1_T_S))).*delta_t;
    
    T_N(i+1)=T_N(i)+0.5.*(k1_T_N+k2_T_N); 
    T_S(i+1)=T_S(i)+0.5.*(k1_T_S+k2_T_S);
    T_E(i+1)=T_E(i)+0.5.*(k1_T_E+k2_T_E);
    
end

%III. Plot it all 
figure(1)
plot(t,T_N,'r')
hold on
plot(t,T_S,'b')
plot(t,T_E,'g')
legend({'T_N','T_E','T_S'})
title('Average Temperature Across Latitudes & Time')
xlabel('time (yrs)')
ylabel('Temperature (K)')
hold off
end
clear all
%% 1.c.
for i=1:1
%I. Setup
delta_t=.1; %yrs 
T_N(1)=260; %Kelvin
T_E(1)=260; %Kelvin
T_S(1)=260; %Kelvin
a_E=.5; %Equatorial Albedo
a_S=.7; %Southern Albedo
e=10E-7; %epsilon value 
S=1250; %flux from sun
D=.5; %factor for 
t=0:delta_t:1000; %generating time steps

a_N=ones(length(t),1).*.7;
for k=501:10000
a_N(k)=a_N(k-1)-.0005; 
if a_N(k) < 0
    a_N(k)=0;
end
end

for i=1:length(t)
    
    a_Nx=a_N(i);
    
    k1_T_N=((1-a_Nx).*S-e.*T_N(i).^4+D.*(T_E(i)-T_N(i))).*delta_t;
    k1_T_S=((1-a_S).*S-e.*T_S(i).^4+D.*(T_E(i)-T_S(i))).*delta_t;
    k1_T_E=((1-a_E).*S-e.*T_E(i).^4+D.*((2.*T_E(i))-T_N(i)-T_S(i))).*delta_t;
    
    k2_T_N=((1-a_Nx).*S-e.*(T_N(i)+k1_T_N)^4+D.*((T_E(i)+k1_T_E)-(T_N(i)-k1_T_N))).*delta_t;
    k2_T_S=((1-a_S).*S-e.*(T_S(i)+k1_T_S)^4+D.*((T_E(i)+k1_T_E)-(T_S(i)-k1_T_S))).*delta_t;
    k2_T_E=((1-a_E).*S-e.*(T_E(i)+k1_T_E)^4+D.*((2.*(T_E(i)+k1_T_E))-(T_N(i)-k1_T_N)-(T_S(i)-k1_T_S))).*delta_t;
    
    T_N(i+1)=T_N(i)+0.5.*(k1_T_N+k2_T_N); 
    T_S(i+1)=T_S(i)+0.5.*(k1_T_S+k2_T_S);
    T_E(i+1)=T_E(i)+0.5.*(k1_T_E+k2_T_E);
   
end

%III. Plot it all 
figure (2)
plot(t,T_N(1:10001),'r')
hold on
plot(t,T_S((1:10001)),'b')
plot(t,T_E(1:10001),'g')
legend({'T_N','T_E','T_S'})
title('Average Temperature Across Latitudes & Time')
xlabel('time (yrs)')
ylabel('Temperature (K)')
hold off

disp('Glacial Melting at the North Pole could explain this phenomenon.')
end
%clear all
%% 1.d.
for i=1:1
 %I. Setup
delta_t=.1; %yrs 
T_N(1)=260; %Kelvin
T_E(1)=260; %Kelvin
T_S(1)=260; %Kelvin
a_N=.7; %Northern Albedo
a_E=.5; %Equatorial Albedo
a_S=.7; %Southern Albedo
e=10E-7; %epsilon value 
D=.5; %factor for 
t=0:delta_t:1000; %generating time steps

S=1250; %flux from sun
flux_t=1:1000;
dv=.0005.*S.*(1/1000);
flux=sin(flux_t).*dv;

%ok now all i need to do it get flux in sets of 10 for each value
%and then we're good 


for i=1:(length(t)-1)

    k1_T_N=((1-a_N).*S-e.*T_N(i).^4+D.*(T_E(i)-T_N(i))).*delta_t;
    k1_T_S=((1-a_S).*S-e.*T_S(i).^4+D.*(T_E(i)-T_S(i))).*delta_t;
    k1_T_E=((1-a_E).*S-e.*T_E(i).^4+D.*((2.*T_E(i))-T_N(i)-T_S(i))).*delta_t;
    
    k2_T_N=((1-a_N).*S-e.*(T_N(i)+k1_T_N)^4+D.*((T_E(i)+k1_T_E)-(T_N(i)-k1_T_N))).*delta_t;
    k2_T_S=((1-a_S).*S-e.*(T_S(i)+k1_T_S)^4+D.*((T_E(i)+k1_T_E)-(T_S(i)-k1_T_S))).*delta_t;
    k2_T_E=((1-a_E).*S-e.*(T_E(i)+k1_T_E)^4+D.*((2.*(T_E(i)+k1_T_E))-(T_N(i)-k1_T_N)-(T_S(i)-k1_T_S))).*delta_t;
    
    T_N(i+1)=T_N(i)+0.5.*(k1_T_N+k2_T_N); 
    T_S(i+1)=T_S(i)+0.5.*(k1_T_S+k2_T_S);
    T_E(i+1)=T_E(i)+0.5.*(k1_T_E+k2_T_E);
    
end

%III. Plot it all 
figure(1)
plot(t,T_N,'r')
hold on
plot(t,T_S,'b')
plot(t,T_E,'g')
legend({'T_N','T_E','T_S'})
title('Average Temperature Across Latitudes & Time')
xlabel('time (yrs)')
ylabel('Temperature (K)')
hold off
end   
clear all
%% 1.e. 
for i=1:1
disp('Cloud reflectivity and its effect on albedo is a big issue in climate')
disp('models due to the range of ways that clouds can either absorb more')
disp('sunlight or reflect it back out. It is thought that on average, more')
disp('clouds cause a high albedo, which would probably cause temperatures to')
disp('decrease, since the light is being reflected back out. Another thing that')
disp('is missing from this model is the latitudinal change in albedo (other than')
disp('North, South, and Equatorial latitudes). There is a much more refined')
disp('resolution of albedo that varies even from gridbox to gridbox in ESMs that')
disp('would change the outcome of this toy model (though I am not sure in what')
disp('direction. Another thing that could potentially change the outcome of this')
disp('model that was not included is the actual temperature difference between')
disp('latitudes- its obvious that there would be a drastic difference in')
disp('starting temperature across the different latitudes, which could then be')
disp('even better refined with grid boxes for subtropical and subpolar')
disp('latitudes.')
end