% minimize wing spar weight subject to stress constraints at manuever
clear all;
close all;

% carbon fiber values from http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp
Nelem = 30;
L = 7.5; % semi-span in meters
rho = 1600; % density of standard carbon fiber, kg/m^3
yield = 600e6; % tensile strength of standard carbon fiber, Pa
E = 70e9; % Young's modulus, Pa
W = 0.5*500*9.8; % half of the operational weight, N
%f = force(W,L,Nelem,Xi); % loading at manueuver

% define function and constraints
fun = @(r) SparWeight(r, L, rho, Nelem);
nonlcon = @(r) WingConstraints(r, L, E, yield, Nelem, W);
lb = 0.01*ones(2*(Nelem+1),1);
up = 0.05*ones(2*(Nelem+1),1);
A = zeros(Nelem+1,2*(Nelem+1));
b = -0.0025*ones(Nelem+1,1);
for k = 1:(Nelem+1)
    A(k,k) = 1.0;
    A(k,Nelem+1+k) = -1.0;
end

% define initial guess (the nominal spar)
r0 = zeros(2*(Nelem+1),1);
r0(1:Nelem+1) = 0.04625*ones(Nelem+1,1);
r0(Nelem+2:2*(Nelem+1)) = 0.05*ones(Nelem+1,1);

options = optimset('GradObj','on','GradConstr','on', 'TolCon', 1e-4, ...
    'TolX', 1e-8, 'Display','iter', 'Algorithm', 'SQP'); %, 'DerivativeCheck','on');
[ropt,fval,exitflag,output] = fmincon(fun, r0, A, b, [], [], lb, up, ...
    nonlcon, options);

% plot optimal radii
r_in = ropt(1:Nelem+1);
r_out = ropt(Nelem+2:2*(Nelem+1));
x = [0:L/Nelem:L].';
figure
plot(x, r_in, '-ks');
hold on;
plot(x, r_out, '--ks');
xlabel('Position on Wing Spar (m)')
ylabel('Distance from the Center of the Spar')
title("Radii of the Spar Along the Wing for the Optimal Design")
legend("Inner Radius", "Outer Radius")

% display weight and stress constraints
[f,~] = fun(ropt)
[c,~,~,~] = nonlcon(ropt)

%%
Iyy = CalcSecondMomentAnnulus(r_in, r_out);

[sigma_mean, Stan_dev] = Rand_Var_Def(r_in, r_out, L, E, Iyy, Nelem, W);

figure
hold on
plot(x,sigma_mean,"-")
plot(x,sigma_mean+6*Stan_dev,"--o")
plot(x,sigma_mean-6*Stan_dev,"--*")
xlabel('Position on Wing Spar (m)')
ylabel('Stess (Pa)')
title("Stress Along the Wing for the Optimal Design")
legend("Mean Stress", "Mean + 6 Std", "Mean - 6 Std")