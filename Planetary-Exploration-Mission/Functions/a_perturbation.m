function [a_tot]=a_perturbation(t,x)

% AUTHOR: 
 % Lorenzo Dionigi

% SUPERVISOR:
 % Prof. Camilla Colombo

sma=x(1); %sma is semi major axis
e=x(2);
i=x(3);
Omega=x(4);
omega=x(5);
theta=x(6);
J2=astroConstants(9);
Re = astroConstants(23);
mu = astroConstants(13);
p=sma*(1-e^2);
r=p/(1+e*cos(theta));


P = 4.5*1e-6;
r_s=astroConstants(2); %AU
C_r=1.0;
A_M=5; %A/M = [m^2/kg]


J2 = astroConstants(9);
SRP = P*r_s^2/((abs(r_s))^2)* C_r*A_M; %[N/kg]
SRP= SRP*1e-3; 
r_e_s =[1; 0;0]*r_s;%radious earth - sun
r_sc_s =r-r_e_s;%radious s/c - sun
a_SRP = -SRP*r_sc_s/norm(r_sc_s);

a_SRP= a_SRP;



a_J2 = -(3*J2*mu*Re^2)/(2*r^4)*[1-3*(sin(i))^2*(sin(theta+omega))^2;
                                                       (sin(i))^2*sin(2*(theta+omega));
                                                       sin(2*i)*sin(theta+omega)];
a_tot=a_J2+a_SRP;
end

