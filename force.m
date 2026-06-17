function f = force(W,L,Nelem,Xi)
% Computes the loading on a wing for the wing-spar problem in Newtons
% Inputs:
%   W - half of the operational weight, N
%   L - length of the beam
%   Nelem - number of finite elements to use
%   Xi - Random variables that change the loading 
% Outputs:
%   f - force per unit length along the beam axis x
%--------------------------------------------------------------------------
fnom = (2*(2.5*W)/(L^2))*(L:-L/Nelem:0).';

Delta_f = zeros(size(fnom));

for n = 1:length(Xi)
    Delta_f = Delta_f + Xi(n)*cos(((2*n-1)*pi*(0:L/Nelem:L).')/(2*L));
end

f = fnom + Delta_f;

end