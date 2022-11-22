

function  MetLV() ;

%Metropolis MCMC for fitting parameters in the Lotka-Volterra (LV) model

close all;

%% LV model parameters
T = 50; %Time period for LV model

P=[0.2,0.05,0.5]; %parameter vector for k1,k2,k3

y0=[5;10]; %intial values

tspan = [0:2:T]; %solutions computed at these timepoints

%% Metropolis parameters

N = 2e4;  %number of steps in the Metropolis algorithms

burntime = N - 1e4; % we only save results after burntime steps

sigma = 1e-2*[3,1,1]; %sigma for computing new proposed values

sigma_chi = 4.; %sigma occuring in the error function chi^2


%% simulate and plot data by running the LV model

[t,data] = ode45(@dy_o_dt,tspan,y0,[],P);


plot(t, data(:,1),'r', t,data(:,2),'g','LineWidth',3)
legend('X','Y');
xlabel('time');
ylabel('concentration');

%% Metropolis parameter fitting

theta_current = 2e-1*ones(1,3); %first guess for theta

[t,y] = ode45(@dy_o_dt,tspan,y0,[],theta_current);

chi_current = sum((data(:) - y(:)).^2); %sigma_chi will be introduced later

% allocate space for saving theta, chi
theta=zeros(N-burntime,3);
chi=zeros(N-burntime,1);

% main iteration
for n = 1:N
    
    
    theta_proposed = theta_current+ sigma.* randn(1,3);
    
    [t,y] = ode45(@dy_o_dt,tspan,y0,[],theta_proposed); %run model for proposed theta
       
    chi_proposed = sum((data(:) - y(:)).^2);
    
    
    ratio=exp((-chi_proposed+chi_current)./(2*sigma_chi^2));
    
    
    if any(theta_proposed<0) %ensure parameters are positive
        ratio=0;
    end
    
    
    if rand < ratio
        theta_current=theta_proposed;
        chi_current=chi_proposed;
    end
    
    
    if n > burntime %save theta and chi after burntime
        theta(n-burntime,:)=theta_current;
        chi(n-burntime)=chi_current;
    end
    
end % end of main iteration

%% Output plots

% plot theta for successive iterations (>burntime>
figure (2)
plot(1:N-burntime,theta);
legend('P(1)','P(2)','P(3)')

% Determine theta for which chi^2 is minimal
minchi=min(chi)
minchi_index=find(chi==minchi,1);
theta(minchi_index,:)

% plot histograms of the parameter distributions
figure (3)

subplot(1,3,1)
hist(theta(:,1));
subplot(1,3,2)
hist(theta(:,2));
subplot(1,3,3)
hist(theta(:,3));

titlestring=sprintf('P=%.2f , %.2f , %.2f',P(1),P(2),P(3))
sgtitle(titlestring);

% plot distribution of error for the first two parameters
figure (4)

x=theta(:,1);
y=theta(:,2);
z=chi;

xv = linspace(min(x), max(x), 40);
yv = linspace(min(y), max(y), 40);
[X,Y] = meshgrid(xv, yv);

Z = griddata(x,y,z,X,Y);
surf(X, Y, Z);
grid on

% plot model run with the fitted parameters
figure (5)
[t,data]=ode45(@dy_o_dt,tspan,y0,[],theta(minchi_index,:));


plot(t, data(:,1),'r', t,data(:,2),'g','LineWidth',3)
legend('X','Y');
xlabel('time');
ylabel('concentration');
title('model with fitted parameters');




end


%% LV model equations
function dy = dy_o_dt(t,y,P)

dy=[ P(1)*y(1) - P(2)*y(1)*y(2);
    P(2)*y(1)*y(2) - P(3)*y(2)];

end





