function LVfit_leastsquares

%parameters of the Lotka-Volterra model
P = [0.2,0.05,0.5];

%initial paramter guesses for fitting the model to simulated data
Pest = [0.25,0.1,1];

%initial conditions
y0 = [5;10];

%time span
tdata = [0:2:50];

%we simulate the model
[t,data] = ode45(@dy_o_dt,tdata,y0,[],P);

%and add some noise with standard deviation sd
sd = 1; 
data = data+sd*randn(size(data));

%fit the simulated data to the original model
Pest = lsqnonlin(@objectivefunction,Pest,[0,0,0],[],[],t,data,y0);

%compute the estimated solution using the fitted Parameters
[t,est] = ode45(@dy_o_dt,tdata,y0,[],Pest);

%plot results

figure (1)

p=plot(t,data(:,1),'-r',t,data(:,2),'-b',t,est(:,1),'--r',t,est(:,2),'--b');

set(p,'linewidth',3); 
legend('y1','y2','y1 fit','y2 fit');
title('Lotka Volterra fit using lsqnonlin');


%% Lotka Volterra model
    function dy = dy_o_dt(t,y,P)
        
        dy=[ P(1)*y(1) - P(2)*y(1)*y(2);
             P(2)*y(1)*y(2) - P(3)*y(2)];
            
    end

%% objective function which is to be minimised
    function e = objectivefunction(Pest,t,data,y0)
        
        [t,y]=ode45(@dy_o_dt,t,y0,[],Pest);
   
        e = y - data; %error, note that lsqnonlin does not expect the squared error
    
    end



end

