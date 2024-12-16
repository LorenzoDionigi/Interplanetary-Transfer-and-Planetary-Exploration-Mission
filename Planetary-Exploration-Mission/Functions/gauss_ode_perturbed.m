function dy = gauss_ode_perturbed(t, y, perturbation_ID, parameters)
% Builds ODE system in Keplerian Parameter
%
% INPUT
% fn            ode system                          [-]
%
% OUTPUT
% t:            time                                [s]
% y:            state                               [km, km/s]

% AUTHOR: 
 % Lorenzo Dionigi

% SUPERVISOR:
 % Prof. Camilla Colombo


dy = zeros(size(y));
a = y(1);
e = y(2);
i = y(3);
% OM = y(4); Not needed
om = y(5);
th = y(6);

% Constants
mu = parameters.mu;
p = a*(1-e^2);
r = p/(1 + e*cos(th));
b = a*sqrt(1-e^2);
n = sqrt(mu/a^3);
h = n*a*b;

% Accelerations
a_tot = compute_disturb(t, y, perturbation_ID, parameters);
a_r = a_tot(1);
a_s = a_tot(2);
a_w = a_tot(3);
% Update
dy(1) = 2*a^2/h * (e*sin(th)*a_r + p/r*a_s);
dy(2) = 1/h*(p*sin(th)*a_r + ((p + r)*cos(th) + r*e)*a_s);
dy(3) = r*cos(th + om)*a_w/h;
dy(4) = r*sin(th+om)*a_w/(h*sin(i));
dy(5) = 1/(h*e)*(-p*cos(th)*a_r + (p + r)*sin(th)*a_s) - (r*sin(th+om)*cos(i)*a_w)/(h*sin(i));
dy(6) = h/r^2 + 1/(e*h)*(p*cos(th)*a_r - (p+r)*sin(th)*a_s);

end




