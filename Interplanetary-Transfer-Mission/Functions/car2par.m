function [a,e,i,OM,om, th ] = car2par(rr, vv,mu)

% [a,e,i,OM,om, th ] = car2par(rr, vv,mu)
% From cartesian to parametric parameters of an orbit.

% INPUT:
 % rr               position vector [Km]
 % vv               velocity vector [Km/s]
 % mu               planetari gravitational constant [Km^3 / s^2]

% OUTPUT:
 % a                major semiaxis [Km]
 % e                eccentricity [-]
 % i                inclination [rad]
 % OM               RAAN [rad]
 % om               pericenter anomaly [rad]
 % theta            thrue anomaly [rad]

% CONTRIBUTORS: 
 % Lorenzo Dionigi, Bouchra Bouras, Giuseppe Antonio Zito, Francesco Paolo Vacca

% SUPERVISOR:
 % Prof. Camilla Colombo

r_mod=norm(rr);
v_mod=norm(vv);

a=1/(2/r_mod - (v_mod^2 / mu));

h= cross(rr,vv);
h_mod=norm(h);
 
e_vect= cross(vv,h)/mu - rr/r_mod;
e=norm(e_vect);

i=acos(h(3)/h_mod);

n=cross([0;0;1],h)/norm(cross([0;0;1],h));

if(n(2)>=0)
   OM=acos(n(1));
else
   OM=2*pi- acos(n(1));
end

if(dot(e_vect,[0;0;1])>=0)
   om=acos(dot(n,e_vect)/e);
else
   om=2*pi- acos(dot(n,e_vect)/e);
end

if(dot(rr,vv)>=0)
   th=acos(dot(rr,e_vect)/(r_mod*e));
else
   th=2*pi- acos(dot(rr,e_vect)/(r_mod*e));
end
