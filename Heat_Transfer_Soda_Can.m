close all;

% Define parameters
h = 6.8; % Convective heat transfer coefficient
r = 65.96 / 1000; % radius
height = 103.3 / 1000; % height  
A_s = 2 * pi * r * height + pi * r^2; % Surface area
rho = 999; % Density
c_p = 4.184 * 1000; % Specific heat [J/kg*K]
k = 22.7e-3; % Air
k_water = 598e-3; % Water 
V = pi * r^2 * height; % Volume

T_init = 25; % Initial temperature
T_infinity = -15; % Ambient temperature
theta_i = T_init - T_infinity; % Initial temperature difference

Tau = (1 / (h * A_s)) * (rho * c_p * V);

% Define time vector
t = linspace(0, 5 * 60 * 60, 10000);

% Calculate Ts(t)
T_s = theta_i * exp(-t / Tau) + T_infinity;
q_s = h * (T_s - T_infinity);

figure;
plot(t, T_s, 'b', 'LineWidth', 2);
xlabel('Time');
ylabel('Temperature (Ts)');
title('Temperature Variation');
grid on;

% Create the geometry
gdm = [
    3   4   0   0   r   r  -height/2  height/2  height/2  -height/2;
];

g = decsg(gdm');

model = createpde("thermal","transient");
geometryFromEdges(model, g);

figure;
pdegplot(model,'EdgeLabels','on');
ylim([-0.06, 0.06]); % Set y-axis limits
axis equal; % Set aspect ratio to be equal

thermalProperties(model,"ThermalConductivity", k_water,...
                        "MassDensity", rho,...
                        "SpecificHeat", c_p);

% Define the time-dependent Neumann boundary condition function
thermalBC(model, 'Edge', 2, 'HeatFlux', @(location,state) mygfun(location,state));
thermalBC(model, 'Edge', 3, 'HeatFlux', @(location,state) mygfun(location,state));

% Initial Condition
thermalIC(model, T_init); % Use T_init instead of 25

% Mesh Generation
generateMesh(model);
figure
pdemesh(model)
title("Mesh with Quadratic Triangular Elements")

% Solve the model
thermalresults = solve(model, t);

% Extract heat flux
[qx,qy] = evaluateHeatFlux(thermalresults);
temperature_values = thermalresults.Temperature;

% Plot the temperature distribution
pdeplot(thermalresults.Mesh.Nodes, thermalresults.Mesh.Elements, ...
        'XYData', temperature_values(:, end), ...
        'Contour', 'on', ...
        'ColorMap', 'hot');
title('Temperature Distribution');

% Add labels and adjust axis as needed
xlabel('x');
ylabel('y');
caxis([-20, 25]);

% Define the function mygfun
function funOut = mygfun(location, state)
    h = 6.8; % Convective heat transfer coefficient
    r = 65.96 / 1000; % radius
    height = 103.3 / 1000; % height  
    A_s = 2 * pi * r * height + pi * r^2; % Surface area
    rho = 999; % Density
    c_p = 4.184 * 1000; % Specific heat [J/kg*K]
    k = 22.7e-3; % Air
    k_water = 598e-3; % Water 
    V = pi * r^2 * height; % Volume
    
    T_init = 25; % Initial temperature
    T_infinity = -15; % Ambient temperature
    theta_i = T_init - T_infinity; % Initial temperature difference
    
    Tau = (1 / (h * A_s)) * (rho * c_p * V);
    time = state.time;
    
    % Calculate Ts(t)
    T_s = theta_i * exp(-time / Tau) + T_infinity;
    q_s = h * (T_s - T_infinity);
    
    funOut = -q_s;
end