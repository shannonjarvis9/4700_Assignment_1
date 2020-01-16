clear all
clearvars
clearvars -GLOBAL
close all
format shorte

% Initialize fundamental constants 
global C
C.temp = 300;                       % Initial temperature 
C.kb = 1.3806504e-23;               % Boltzmann constant
C.m_0 = 9.10938215e-31;             % Electron mass
C.m_e = 0.26*C.m_0;                 % Effective mass 

numAtoms = 1;                       % Number of particles 
delta_t = 0.5;                      % Time step 
numSteps = 1000;                    % Number of time stpes 

% Thermal Velocity 
vth = sqrt((2 * C.kb * C.temp) / C.m_e);

% Initialize the particle position 
x = 200*rand(1,numAtoms); %Size of region in x-direction is 200nm 
y = 100*rand(1,numAtoms); %Size of region in y-direction is 100nm


% Assign each particle fixed velocity given by v_th (found in part 1 using
% temp)
% std0 = sqrt(C.kb * C.temp / C.m_e);
% Vx = std0 * randn(1,numAtoms);
% Vy = std0 * randn(1,numAtoms);
Vx = 0.5; %Change later 
Vy = 0.5;

% Assign a random direction 
phi = rand(1, numAtoms)* 2 * pi; 

%Move electrion in that direction 
x2 = x + cos(phi).*Vx*delta_t;
y2 = y + sin(phi).*Vy*delta_t;

Px = zeros(1, numSteps);
Py = zeros(1, numSteps); 

Px(1,1) = x;
Py(1,1) = y;
Px(1,2) = x2;
Py(1,2) = y2;
%loop over time in steps, plot position of points at each time step 
%Don't loop over particles, just keep the current and prev posiiton 
for k1= 3:numSteps
    x_prev = x2;
    y_prev = y2;
    
    x2 = x_prev + cos(phi).*Vx*delta_t;
    y2 = y_prev + sin(phi).*Vy*delta_t;

    
    %Add boundary conditions 
    if x2> 200 
        x2 = x2 - 200;
    end 
    
    if x2 < 0
        x2 = x2 + 200;
    end 
    
    if y2 > 100
        y2 = y_prev + sin(phi).*Vy*delta_t;
    end 
    
    j = state(:,2) > height;
    state(j,2) = 2*height - state(j,2);
    state(j,4) = -state(j,4);
    
    j = state(:,2) < 0;
    state(j,2) = -state(j,2);
    state(j,4) = -state(j,4);
    
    Px(1,k1) = x2;
    Py(1,k1) = y2;
end 

plot(Px, Py)
% Add rules for boundary conditions 
%updated_x = x + (delta_t)(velocity) 
%Boundary conditions rules 
