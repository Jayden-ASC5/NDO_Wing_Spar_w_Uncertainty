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

%% Compute the mean and standard deviation of stress along the beam
%% Initializing the random variables
r = 4;                        % Number of random variables
Xi = zeros(r,1);              % Random variable array
Stan_dev_f = zeros(r,1);      % Standard deviation array for each random variable

% Define the nominal force distribution (fnom) along the beam
fnom = (2*(2.5*W)/(L^2)) * ((L:-L/Nelem:0).');  % Linear force distribution

% Calculate standard deviations of the forces for each random variable
for i = 1:length(Xi)
    Stan_dev_f(i,1) = fnom(1) / (10 * i);  % Scaling std deviation by the index
end


Stress_root = zeros(6,1);

%% Gauss-Hermite quadrature
xi_full = [0 -0.707107 -1.22474487139 -1.650680 -2.020183 -2.350605; 
           0 0.707107 0 -0.524648 -0.958572 -1.335849; 
           0 0 1.22474487139 0.524648 0 -0.436077;
           0 0 0 1.650680 0.958572 0.436077;
           0 0 0 0 2.020183 1.335849;
           0 0 0 0 0 2.350605]; % Quadrature points (roots of Hermite polynomial)

wts_full = [sqrt(pi) 0.886227 0.295408975151 0.081313 0.019953 0.004530; 
            0 0.886227 1.1816359006 0.804914 0.393619 0.157067; 
            0 0 0.295408975151 0.804914 0.945309 0.724629;
            0 0 0 0.081313 0.393619 0.724629;
            0 0 0 0 0.019953 0.157067;
            0 0 0 0 0 0.004530] ./ sqrt(pi); % Quadrature weights

for a = 1:6
    xi = xi_full(1:a,a);
    wts = wts_full(1:a, a);

    % Initialize mean and squared mean arrays for stress calculation
    mean_f = zeros(Nelem+1,1);    % Array to store mean stress values
    mean_f2 = zeros(Nelem+1,1);   % Array to store squared mean stress values

    %% Computing Gaussian-Hermite Quadrature to find the mean and variance of stress
    for i1 = 1:size(xi,1)
        pt1 = sqrt(2) * Stan_dev_f(1,1) * xi(i1);  % Evaluate random variable 1
        for i2 = 1:size(xi,1)
            pt2 = sqrt(2) * Stan_dev_f(2,1) * xi(i2);  % Evaluate random variable 2
            for i3 = 1:size(xi,1)
                pt3 = sqrt(2) * Stan_dev_f(3,1) * xi(i3);  % Evaluate random variable 3
                for i4 = 1:size(xi,1)
                    pt4 = sqrt(2) * Stan_dev_f(4,1) * xi(i4);  % Evaluate random variable 4
                    
                    % Combine random variables into a vector
                    Xi = [pt1; pt2; pt3; pt4];
    
                    % Compute force distribution along the beam for the given random variables
                    f = force(W, L, Nelem, Xi);
                    
                    % Compute beam displacements under the given forces
                    u = CalcBeamDisplacement(L, E, Iyy, f, Nelem);
                    
                    % Compute stresses based on displacements
                    sigma = CalcBeamStress(L, E, r_out, u, Nelem);
    
                    % Update the mean stress and mean squared stress using quadrature weights
                    mean_f = mean_f + wts(i1) * wts(i2) * wts(i3) * wts(i4) * sigma;
                    mean_f2 = mean_f2 + wts(i1) * wts(i2) * wts(i3) * wts(i4) * sigma.^2;
                end
            end
        end
    end

    % Plot the current value of Stress_root
    % Compute the mean stress
    sigma_mean = mean_f;

    % Compute the standard deviation of stress
    Stan_dev = (mean_f2 - (mean_f).^2) .^ (1/2);

    Stress = sigma_mean + 6 * Stan_dev;
    
    Stress_root(a) = Stress(1);
    
end
x = (0:L/Nelem:L).';
plot(x,f)
xlabel('Position on Wing Spar (m)')
ylabel('Force (N)')
title("Force Along the Wing")

plot(1:6, (Stress_root(1:6)-Stress_root(6)), '-o', 'LineWidth', 1.5);
xlabel('Number of Quadrature points')
ylabel('Stress (Pa)')
title('Comparison Between the Stress at the Root with Different Quadrature Points')

% Calculate percent error between Stress_root(1:5) and Stress_root(6)
percent_error = zeros(5, 1); % Preallocate array for storing percent errors
for i = 1:5
    percent_error(i) = abs((Stress_root(i) - Stress_root(6)) / Stress_root(6)) * 100;
end

% Display the results
disp('Percent Error between Stress_root(1:5) and Stress_root(6):');
for i = 1:5
    fprintf('Iteration %d: %.2f%%\n', i, percent_error(i));
end


% Compute the mean stress
sigma_mean = mean_f;

% % Compute the standard deviation of stress
% Stan_dev = (mean_f2 - (mean_f).^2) .^ (1/2);
% % Plotting the stress distribution along the wing spar
% figure  % Create a new figure for plotting
% hold on % Retain plots for overlaying multiple data series
% plot(x, sigma_mean, "-", 'LineWidth', 1.5)           % Plot mean stress
% plot(x, sigma_mean + 6 * Stan_dev, "--o", 'LineWidth', 1.2) % Mean + 6 std dev
% plot(x, sigma_mean - 6 * Stan_dev, "--*", 'LineWidth', 1.2) % Mean - 6 std dev
% xlabel('Position on Wing Spar (m)')                  % Label for x-axis
% ylabel('Stress (Pa)')                                % Label for y-axis
% title("Stress Along the Wing for the Nominal Design") % Plot title
% legend("Mean Stress", "Mean + 6 Std", "Mean - 6 Std") % Add legend for the plots
