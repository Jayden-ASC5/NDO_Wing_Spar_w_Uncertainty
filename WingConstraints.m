function [c, ceq, dcdx, dceqdx] = WingConstraints(x, L, E, yield, Nelem, W)
% Computes the nonlinear inequality constraints for the wing-spar problem
% Inputs:
%   x - the DVs; x(1:Nelem+1) inner and x(Nelem+2:2*(Nelem+1) outer radius
%   L - length of the beam
%   E - longitudinal elastic modulus
%   force - force per unit length along the beam axis x
%   yield - the yield stress for the material
%   Nelem - number of finite elements to use
% Outputs:
%   c, ceq - inequality (stress) and equality (empty) constraints
%   dcdx, dceqdx - Jacobians of c and ceq
%--------------------------------------------------------------------------
% assert( size(force,1) == (Nelem+1) );
assert( size(x,1) == (2*(Nelem+1)) );

c = CalcInequality(x);
ceq = [];
dcdx = zeros(2*(Nelem+1),Nelem+1);
dceqdx = [];
for k = 1:2*(Nelem+1)
    xc = x;
    xc(k) = xc(k) + complex(0.0, 1e-30);
    dcdx(k,:) = imag(CalcInequality(xc))/1e-30;
end 

    function cineq = CalcInequality(x)
        % compute the displacements and the stresses
        r_in = x(1:Nelem+1);
        r_out = x(Nelem+2:2*(Nelem+1));
        Iyy = CalcSecondMomentAnnulus(r_in, r_out);

        [sigma_mean, Stan_dev] = Rand_Var_Def(r_in, r_out, L, E, Iyy, Nelem, W);

        cineq = (sigma_mean + 6.*(sqrt(Stan_dev)))./yield - ones(Nelem+1,1);

    end
end

