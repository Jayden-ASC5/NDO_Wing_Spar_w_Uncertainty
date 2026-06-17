function [xi, wts] = GaussHermiteQuadrature(order)
    % Computes Gaussian-Hermite quadrature points and weights
    % Inputs:
    %   order: Number of quadrature points
    % Outputs:
    %   xi: Quadrature points (roots of Hermite polynomial)
    %   wts: Quadrature weights

    % Validate the order
    if order < 1
        error('Order must be at least 1.');
    end

    % Generate the Hermite polynomial coefficients using recurrence
    H = zeros(order + 1, order + 1);  % Matrix to store coefficients
    H(1, end) = 1;  % H_0(x) = 1
    if order > 1
        H(2, end-1:end) = [2, 0];  % H_1(x) = 2x
    end
    for n = 2:order
        % H_{n+1}(x) = 2x * H_n(x) - 2n * H_{n-1}(x)
        H(n + 1, :) = 2 * [H(n, 2:end), 0] - 2 * (n - 1) * H(n - 1, :);
    end

    % Roots of the Hermite polynomial give quadrature points
    xi = sort(roots(H(order + 1, :)));

    % Compute weights
    wts = zeros(order, 1);
    for i = 1:order-1
        % Evaluate the derivative of H_n(x) at each root
        H_prime = polyder(H(order + 1, :));
        wts(i) = 2^(order - 1) * factorial(order) * sqrt(pi) / ...
                 (polyval(H_prime, xi(i))^2);
    end
end
