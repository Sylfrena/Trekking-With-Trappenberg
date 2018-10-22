% ode for wilson equations 

%% parameters of the model: 1=K,R; 2=Ca,T; 3=KCa,H; 4=Na; 

function ydot=wilson_ode(t,y,flag,I_ext) 
    g(1)=26;
    g(2)=2.25; 
    g(3)=9.5; 
    g(4)=1; 
    g=g’; 
    
    E(1)=-.95; 
    E(2)=1.20; 
    E(3)=E(1); 
    E(4)=.50; 
    E=E’; 

    tau(1)=1/4.2; 
    tau(2)=1/14; 
    tau(3)=1/45; 
    tau=tau’; 

    V=y(4);

    y3=y(1:3); 

    x0(1)=1.24 + 3.7*V + 3.2*V^2;
    x0(2)=4.205 + 11.6*V + 8 *V^2;
    x0(3)=3*y(2);
    x0=x0’; 

    ydot=tau.*(x0-y3); 

    y3(4)=17.8 + 47.6*V +33.8*V^2; 

    I=g.*y3.*(y(4)-E); 
    
    ydot(4)=I_ext-sum(I); 
return