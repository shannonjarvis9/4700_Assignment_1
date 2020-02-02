clear all 
clear

% Initialize fundamental constants 
global C
C.temp = 300;                       % Initial temperature 
C.kb = 1.3806504e-23;               % Boltzmann constant
C.m_0 = 9.10938215e-31;             % Electron mass
C.m_e = 0.26*C.m_0;                 % Effective mass 

x_max = 200e-9;                     %maximum x dimension
y_max = 100e-9;                     %maximum y dimension

delta_t = .01e-12;                  % Time step 
numSteps = 100;                     % Number of time stpes 
numAtoms = 1000;                       % Number of particles 
numPlot = 20;                       % Number of particles to plot 

T_avg = zeros(numSteps,1);
time = linspace(0,delta_t*numSteps, numSteps);


% Thermal Velocity = 1.870192676075498e+05
v_th = sqrt((2 * C.kb * C.temp) / C.m_e);

%mean free path = 3.740385352150996e-08
tau = 0.2e-12;
mfp= v_th*tau;


% Initialize the particle position
x = x_max*rand(numAtoms,1);     %Size of region in x-direction is 200nm 
y = y_max*rand(numAtoms,1);     %Size of region in y-direction is 100nm



%Assign a velocity from the Maxwell Boltzmann Distribution 
Vx = v_th.*rand(numAtoms,1);
Vy = v_th.*rand(numAtoms,1);

figure(3);
hist(sqrt(Vx.^2 + Vy.^2));
title('Distribution of electron velocity');
xlabel('Velocity (m/s)');
ylabel('Frequency');

color = hsv(numPlot);

prob_scatter = 1- exp(-delta_t/0.2e-12);

%loop over time in steps, plot position of points at each time step 
%Don't loop over particles, just keep the current and prev posiiton 
for i = 1:numSteps

    %Rethermalize   
    for j=1:length(x)
        if prob_scatter > rand()
                Vx(j)= v_th.*rand(1,1);
                Vy(j)= v_th.*rand(1,1);
        end
    end
    

    y_prev = y;
    x_prev = x; 
    
    
    %Move electron
    x = x + Vx*delta_t;
    y = y + Vy*delta_t;
    
    %Add boundary conditions
    above_x_bounds = logical(x>=x_max);
    below_x_bounds = logical(x<=0);
    
    above_y_bounds = logical(y>=y_max);
    below_y_bounds = logical(y<=0);
    
    x(above_x_bounds) = x(above_x_bounds) - x_max;
    x_prev(above_x_bounds) = 0;
    
    x(below_x_bounds) = x(below_x_bounds) + x_max;
    x_prev(below_x_bounds) = x_max;

    y(above_y_bounds) = -y(above_y_bounds) + 2*y_max;
    Vy(above_y_bounds) = -Vy(above_y_bounds);
    

    Vy(below_y_bounds) = -Vy(below_y_bounds);
    
    figure(1)
    xlabel('(m)');
    ylabel('(m)');
    title('Plot of the particle trajectory');
    axis ([0 x_max 0 y_max]);
    pause(0.0001)
    hold on;
    
    for j=1:numPlot 
        plot([x_prev(j)';x(j)'], [y_prev(j)';y(j)'], 'color', color(j,:)); 
    end 
    
    Vavg = mean(Vx.^2 + Vy.^2); %it is already squared 
    %v_th^2 = (2 * C.kb * C.temp) / C.m_e);
    Tavg(i) = ( Vavg*C.m_e)/(2*C.kb);
    
    mfp = tau*Vavg;
    %average time between colisions 
    tau_calc = (mfp*numAtoms)/Vavg;

    
end


figure(2)
plot(time,Tavg);
str = sprintf('Plot of the average temperature (average is %d K)', mean(Tavg));
title(str);
xlabel('Time');
ylabel('Temperature (K)');

















































































