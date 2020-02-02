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

delta_t = .02e-12;                  % Time step 
numSteps = 100;                     % Number of time stpes 
numAtoms = 500;                       % Number of particles 
numPlot = 20;                       % Number of particles to plot 

specular = 1; %BC for bottleneck boxes

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


%Check if pqarticle is in the box/not allowed region 
for i=1:numAtoms
    while (y(i)<= 40e-9 || y(i)>= 60e-9) && (x(i)>= 80e-9 && x(i) <= 120e-9)
        %Particle is in the box, generate new position 
        x(i) = x_max*rand(); 
        y(i) = y_max*rand(); 
    end
end 


% Assign a random direction 
%Assign a velocity from the Maxwell Boltzmann Distribution 
Vx = v_th.*rand(numAtoms,1);
Vy = v_th.*rand(numAtoms,1);



color = hsv(numPlot);

prob_scatter = 1- exp(-delta_t/0.2e-12);

%loop over time in steps, plot position of points at each time step 
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
hold on;
rectangle('position', [80e-9 0 40e-9 40e-9]);
rectangle('position', [80e-9 60e-9 40e-9 40e-9]);
hold on;
xlim(axes1,[0 x_max]);
ylim(axes1,[0 y_max]);
xlabel('nm');
ylabel('nm');

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
    
    %Add bottleneck conditions 
     if specular == 1 
         for j=1:numAtoms
            if (y(j)<= 40e-9 || y(j) >= 60e-9) && (x(j)>= 80e-9 && x(j) <= 120e-9)
                Vx(j) = - Vx(j);
                x(j) = x_prev(j);
                y(j) = y_prev(j);
            end
            if (y(j) <= 40e-9 && y(j) >= 60e-9) && (x(j) >= 80e-9 && x(j) <= 120e-9)
                Vy(j) = - Vy(j);
                x(j) = x_prev(j);
                y(j) = y_prev(j);
            end
         end
     else  %diffuse/rethermalize 
            for j=1:numAtoms
                while (y(j)<= 40e-9 || y(j)>= 60e-9) && (x(j)>= 80e-9 && x(j) <= 120e-9)
                    %Particle is in the box, generate new position 
                    x(j) = x_max*rand(); 
                    y(j) = y_max*rand();
                    x(j) = x_prev(j);
                    y(j) = y_prev(j);
                end
            end    
     end    
             
             
 
    figure(1)   
    xlabel('(m)');
    ylabel('(m)');
    title('Plot of the particle trajectory');
    axis ([0 x_max 0 y_max]);
    hold on;
    pause(0.0001)
        
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

%Temperature plot 
figure(2)
plot(time,Tavg);
str = sprintf('Plot of the average temperature (average is %d K)', mean(Tavg));
title(str);


% Density Map
figure(2);
hist3([x y], 'CdataMode', 'auto');
colorbar;
grid on;
title('Final particle position density map')
xlabel('x position (m)');
ylabel('y position (m)');
zlabel('Particle count');
hold off;

% Temperature Map
figure(3);
Vavg_final = Vx.^2 + Vy.^2;
Tavg_final = (Vavg_final.*C.m_e)./(2*C.kb);
scatter3(x, y, Tavg_final);
title('Final temperature distribution map')
xlabel('x position (m)');
ylabel('y position (m)');
zlabel('Temperature (K)');
grid on;

