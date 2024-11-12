function dy = ode_2bp_perturbed( t, y, perturbation_ID, parameters)
%ode_2bp ODE system for the two-body problem (Keplerian motion)
%
% PROTOTYPE
% dy = ode_2bp( t, y, mu )
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
%
% OUTPUT:
% dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position and velocity
r = y(1:3);
v = y(4:6);
% Distance from the primary
rnorm = norm(r);
% Computation of the acceleration's perturbation
a_disturb = compute_disturb(t, y, perturbation_ID, parameters);
% Set the derivatives of the state
dy = [ v
         (-parameters.mu/rnorm^3)*r+a_disturb];
end