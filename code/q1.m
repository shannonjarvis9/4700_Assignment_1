clear
close all

% Initialize fundamental constants 
global C
C.temp = 300;                       % Initial temperature 
C.kb = 1.3806504e-23;               % Boltzmann constant
C.m_0 = 9.10938215e-31;             % Electron mass
C.m_e = 0.26*C.m_0;                 % Effective mass 

x_max = 200e-9;                     %maximum x dimension
y_max = 100e-9;                     %maximum y dimension

delta_t = .02e-12;                  % Time step 
numSteps = 100;                     % Number of time stpes 
numAtoms = 20;                      % Number of particles 



% Thermal Velocity = 1.870192676075498e+05
v_th = sqrt((2 * C.kb * C.temp) / C.m_e);

%mean free path = 3.740385352150996e-08
mfp= v_th*0.2e-12;


% Initialize the particle position
x = x_max*rand(numAtoms,1);     %Size of region in x-direction is 200nm 
y = y_max*rand(numAtoms,1);     %Size of region in y-direction is 100nm


% Assign a random direction 
phi = rand(numAtoms,1)* 2 * pi; 

% Assign each particle fixed velocity given by v_th (found in part 1 using
% temp)
Vx = v_th*cos(phi);
Vy = v_th*sin(phi);

colorarray = rand(numAtoms,1);
T_avg = zeros(numSteps,1);
time = 1:numSteps;

color = hsv(numAtoms);

for i = 1:numSteps
    
    x_prev = x;
    y_prev = y;
    
    %Update the Positions
    x = x + Vx*delta_t;
    y = y + Vy*delta_t;

    %Add boundary conditions
    above_x_bounds = x>=x_max;
    x(above_x_bounds) = x(above_x_bounds) - x_max;
    x_prev(above_x_bounds) = 0;
    
    below_x_bounds = x<=0;
    x(below_x_bounds) = x(below_x_bounds) + x_max;
    x_prev(below_x_bounds) = x_max;
    
    above_y_bounds =  y>=y_max;
    Vy(above_y_bounds) = -Vy(above_y_bounds);
    y(above_y_bounds) = y_max - (y(above_y_bounds)-y_max);
    
    below_y_bounds = y<=0; 
    Vy(below_y_bounds) = -Vy(below_y_bounds);
    
    figure(1)
    xlabel('(m)');
    ylabel('(m)');
    str = sprintf('Plot of the particle trajectory for %d particles', numAtoms);
    title(str);
    axis ([0 x_max 0 y_max]);
    pause(0.001)
    hold on;
    
    for j=1:numAtoms 
 
        plot([x_prev(j)';x(j)'], [y_prev(j)';y(j)'], 'color', color(j,:)); 

    end 
    
    
    Vavg = mean(Vx.^2 + Vy.^2); %it is already squared 
    %v_th^2 = (2 * C.kb * C.temp) / C.m_e);
    Tavg(i) = ( Vavg*C.m_e)/(2*C.kb);

    
end

figure(2)
plot(time,Tavg);
str = sprintf('Plot of the average temperature (average is %d K)', Tavg(i));
title(str);
xlabel('Time');
ylabel('Temperature (K)');