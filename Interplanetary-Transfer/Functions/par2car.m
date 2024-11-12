function [rr,vv] = par2car(a,e,i,OM,om, theta,mu)

% [rr,vv] = par2car(a,e,i,Om,om, theta,mu)
% from parametric to cartesian parameters of an orbit

% INPUT:
 % a                major semiaxis [Km]
 % e                eccentricity [-]
 % i                inclination [rad]
 % OM               RAAN [rad]
 % om               pericenter anomaly [rad]
 % theta            thrue anomaly [rad]

% OUTPUT:
 % rr               position vector [Km]
 % vv               velocity vector [Km/s]
 % mu               planetari gravitational constant [Km^3 / s^2]

 % CONTRIBUTORS: 
 % Giuseppe Antonio Zito

p=a*(1-e^2);
r= p /(1+e*cos(theta));
r_pf=r.*[cos(theta);sin(theta);0];
v_pf=sqrt(mu/p)*[-sin(theta);e+cos(theta);0];

R1=@(a)[cos(a) sin(a) 0; -sin(a) cos(a) 0; 0 0 1];
R2=@(b)[1 0 0; 0 cos(b) sin(b); 0 -sin(b) cos(b)];
R3=@(c)[cos(c) sin(c) 0; -sin(c) cos(c) 0; 0 0 1];
R=R1(om)*R2(i)*R3(OM);
rr=R' * r_pf;
vv=R' * v_pf;
