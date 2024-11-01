function [a, T] = find_semimajoraxis(k,m,omega_e, mu)
T_e=2*pi/omega_e;
T= T_e*m/k;
n= 2*pi/T;
a=(mu/(n^2))^(1/3);
a=a;