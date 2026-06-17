% Define parameters for the beam and material properties
Nelem = 15;            % Number of beam elements (discretization of the spar)
L = 7.5;               % Semi-span of the beam (length) in meters
rho = 1600;            % Density of standard carbon fiber, kg/m^3
yield = 600e6;         % Tensile strength of carbon fiber, Pa
E = 70e9;              % Young's modulus (elastic modulus), Pa
W = 0.5 * 500 * 9.8;   % Half of the operational weight, N (0.5 * 500 kg * g)
x = (0:L/Nelem:L).';   % Position along the length of the beam (spatial discretization)

% Initialize design variables for inner and outer radii of the spar
r0 = zeros(2 * (Nelem+1), 1);           % Initialize radii vector
r0(1:Nelem+1) = 0.0415 * ones(Nelem+1, 1);  % Inner radius, constant for all segments
r0(Nelem+2:2*(Nelem+1)) = 0.05 * ones(Nelem+1, 1); % Outer radius, constant for all segments

% Calculate the mass of the spar based on its geometry
mass = SparWeight(r0, L, rho, Nelem);  % Function to compute spar weight
                                       
% Extract inner and outer radii for stress and stiffness calculations
r_in = r0(1:Nelem+1);                  % Inner radii of the spar
r_out = r0(Nelem+2:2*(Nelem+1));       % Outer radii of the spar

% Calculate the second moment of area (Iyy) for the annular cross-section
Iyy = CalcSecondMomentAnnulus(r_in, r_out);  % Moment of inertia for beam elements

% Compute the mean and standard deviation of stress along the beam
[sigma_mean, Stan_dev] = Rand_Var_Def(r_in, r_out, L, E, Iyy, Nelem, W);

% Plotting the stress distribution along the wing spar
figure  % Create a new figure for plotting
hold on % Retain plots for overlaying multiple data series
plot(x, sigma_mean, "-", 'LineWidth', 1.5)           % Plot mean stress
plot(x, sigma_mean + 6 * Stan_dev, "--o", 'LineWidth', 1.2) % Mean + 6 std dev
plot(x, sigma_mean - 6 * Stan_dev, "--*", 'LineWidth', 1.2) % Mean - 6 std dev
xlabel('Position on Wing Spar (m)')                  % Label for x-axis
ylabel('Stress (Pa)')                                % Label for y-axis
title("Stress Along the Wing for the Nominal Design") % Plot title
legend("Mean Stress", "Mean + 6 Std", "Mean - 6 Std") % Add legend for the plots
