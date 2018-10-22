%% Integration of Wilson model with multistep ODE solver 

clear; clf;

%% Integration: 

% Equilibration: no external input; 

y0=zeros(1,4); 
y0(4)=-1; 
param=0; 
I_ext=0; 
tspan=[0 100]; 
[t,y]=ode45(’wilson_ode’,tspan,y0,[],I_ext); 

% Integration with external input; 

y0=y(size(t,1),:); param=0; I_ext=1; tspan=[0 200]; 
[t,y]=ode45(’wilson_ode’,tspan,y0,[],I_ext); 

%% Ploting Results 

plot(t,100*y(:,4)); 
xlabel(’Time’); 
ylabel(’Membrane potential’);