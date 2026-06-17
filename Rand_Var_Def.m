function [sigma_mean, Stan_dev] = Rand_Var_Def(~, r_out, L, E, Iyy, Nelem, W)
% Computes the mean stress and standard deviation of stress
% in a beam using Gaussian-Hermite quadrature for random variables.
%
% Inputs:
%   r_out - Outer radius of the beam cross-section
%   L - Length of the beam
%   E - Young's Modulus (Elastic modulus)
%   Iyy - Second moment of area about y-axis
%   Nelem - Number of beam elements
%   W - Applied load
%
% Outputs:
%   sigma_mean - Mean stress distribution along the beam
%   Stan_dev - Standard deviation of stress distribution along the beam
%--------------------------------------------------------------------------
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

% Initialize mean and squared mean arrays for stress calculation
mean_f = zeros(Nelem+1,1);    % Array to store mean stress values
mean_f2 = zeros(Nelem+1,1);   % Array to store squared mean stress values

%% 3-point Gauss-Hermite quadrature
xi = [-1.22474487139; 0.0; 1.22474487139]; % Quadrature points (roots of Hermite polynomial)
wts = [0.295408975151; 1.1816359006; 0.295408975151] ./ sqrt(pi); % Quadrature weights

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

% Compute the mean stress
sigma_mean = mean_f;

% Compute the standard deviation of stress
Stan_dev = (mean_f2 - (mean_f).^2) .^ (1/2);


end
