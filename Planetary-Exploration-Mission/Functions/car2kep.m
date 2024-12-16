function [a,e,i,OM,om,th] = car2kep(rr,vv, mu)
% Transformation [rr,vv]-->[a,e,i,Ω,ω,θ]
%------------------------------------------------------------
% Input:
% rr = position                                              [km]
% vv = speed                                                 [km/s]
% -----------------------------------------------------------
% Output:
% a = semi-major axis                                       [km]
% e = eccentricity                                          [-]
% i = inclination                                           [rad]
% OM = RAAN                                                 [rad]
% om = pericenter anomaly                                   [rad]
% th0 = initial true anomaly                                [rad]
% thf = final true anomaly                                  [rad]
% dth = true anomaly step size                              [rad]

% AUTHOR: 
 % Lorenzo Dionigi

% SUPERVISOR:
 % Prof. Camilla Colombo

% 1
r = norm(rr);
v = norm(vv);
a = (2/r - v^2/mu)^-1;

% 2
h_vec = cross(rr,vv);
h = norm(h_vec);
e_vec = (cross(vv,h_vec))/mu - rr/r;
e = norm(e_vec);

% 3
i = acos( (h_vec(3)) / h );

% 4 
N = cross([0 0 1],h_vec);
N = N/(norm(N));
if   N(2)>0
    OM = acos( N(1) );
else 
    OM = 2*pi - acos( N(1) );
end

% 5
if e_vec(3) > 0
    om = acos( dot(N,e_vec) / e );
else 
    om = 2*pi - acos( dot(N,e_vec) / e );
end

% 6
V_rad = dot(rr,vv)/r;

if V_rad > 0 
    th = acos( dot(rr,e_vec) / (r * e) );
else 
    th = 2*pi - acos( dot(rr,e_vec) / (r * e) );
end


end


