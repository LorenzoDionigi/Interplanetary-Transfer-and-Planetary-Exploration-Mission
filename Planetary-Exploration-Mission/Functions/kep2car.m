function [r,v]=keptocar(a,e,i,Omega,omega,theta, mu)

% AUTHOR: 
 % Lorenzo Dionigi

% SUPERVISOR:
 % Prof. Camilla Colombo

p=a*(1-e^2);
r_c = p/(1+e*cos(theta));%non è il raggio effettivo ma il raggio della conica
%non sono il raggio e la velocità nel priano cartesiano ma nel piano perifocale
r_vec=[r_c*cos(theta) r_c*sin(theta) 0]'; 
v_vec=sqrt(mu/p)* [-sin(theta) e+cos(theta) 0]';
%trasporto nel piano cartesiano
R1=[cos(Omega)  sin(Omega) 0;
        -sin(Omega) cos(Omega) 0;     %rotazione rispetto asse z
              0                       0           1]; 
R2=[      1                0                 0;
              0            cos(i)       sin(i);   %rotazione rispetto asse x
              0            -sin(i)      cos(i)];      
R3=[cos(omega)  sin(omega) 0;
        -sin(omega) cos(omega) 0;      %rotazione rispetto as
              0                       0          1];

T=R3*R2*R1;
T_inv=T';

%X= T_inv*x --> T*X=x dove X inerziele(cartesiano) e x perifocale
r=T_inv*r_vec;
v=T_inv*v_vec;

