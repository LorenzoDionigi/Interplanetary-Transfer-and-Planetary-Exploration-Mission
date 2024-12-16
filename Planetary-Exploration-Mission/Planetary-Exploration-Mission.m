clc; clearvars; close all;


% This project aims to research the different obstacles when preparing for a planetary explorer
% mission. The main focus is on Earth observation by a spacecraft with a given orbit, and the infuence 
% of orbit perturbations. The perturbations considered are 2 and solar radiation pressure (SRP).
% Characterization of the groundtracks of the spacecraft is used to further analyze the effects of these
% perturbations. The infuence of perturbations is analyzed by propagation through two different
% methods, Cartesian and Gaussian. Lastly a comparison is made between the propagations of the
% given orbit, and the ephemeris of the orbit of COSMOS-2546 with similar orbit.


% CONTRIBUTORS: 
 % Lorenzo Dionigi, Bouchra Bouras, Giuseppe Antonio Zito, Francesco Paolo Vacca

% SUPERVISOR:
 % Prof. Camilla Colombo



addpath('./Functions/');
format long
%% NOMINAL ORBIT
%initial values orbital elements
sma= 1.8886*1e4;%km
e= 0.5795; 
i= deg2rad(64.2795); %rad
OM= pi; %rad
om= pi; %rad
theta0= 0; %rad
x0=[sma, e, i, OM, om, theta0]';


%constants
R_e=astroConstants(23);
mu_E=astroConstants(13);
J2=astroConstants(9);
AU=astroConstants(2);
A_M=5;
Cr=1;
% Compilation of the input struct of constants
parameters.mu = mu_E;
parameters.AU = AU;
parameters.R_e = R_e;
parameters.J2 = J2;
parameters.Cr = Cr;
parameters.A_M = A_M;
parameters.p_Solar = 4.5*1e-6; % [N/m^2]
omega_e=15.04*pi/(180*60*60); %rad/s
T=2*pi*sqrt(sma^3/mu_E); %sec %period of the sc around earth

% n_orbits=24*3600*1/T;
n_orbits=1;

tspan=linspace(0, n_orbits*T, n_orbits*1000); 
% INITIAL ORBIT
[rr,vv]=keptocar(sma, e, i , OM, om, theta0, mu_E);
y0=[rr;vv];
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14);
% Perform the integration
perturbation_ID=0;%0 set unperturbed orbit
odefun= @(t,y) ode_2bp_perturbed(t,y,perturbation_ID,parameters);
[t,Y] = ode113( odefun, tspan, y0, options );

%3D PLOT INIZIAL ORBIT
figure()
Earth_3D()
hold on
plot3(Y(:,1), Y(:,2), Y(:,3),'r','Linewidth',2) 
plot3(Y(1,1), Y(1,2), Y(1,3), 'o', MarkerSize=10, MarkerFaceColor="b") %starting point
xlabel('x [Km]');
ylabel('y [Km]');
zlabel('z [Km]');
legend('','Unperturbed orbit','Initial s/c position')

%% GROUND TRACK INIZIAL ORBIT UNPERTURBED
[~,~] =groundtrack(t,Y,omega_e);

% modifying semi major axis for a repeating ground track
k=10; %number of s/c revolution
m=3; %number of earth revolution
[a, T_GT] = find_semimajoraxis(k,m,omega_e, mu_E);

[r0,v0]=keptocar(a, e, i , OM, om, theta0, mu_E);
y0_repeating_GT=[r0;v0];
% Perform the integration
[t,Y_repeting_GT] = ode113( odefun, tspan, y0_repeating_GT, options );
% REPEATING GROUND TRACK 
[~,~] =groundtrack(t,Y_repeting_GT,omega_e);

%% PERTURBED INITIAL ORBIT 
% perturbations: J2 and SRP
% par:cR= 1.0 A/M = 5.0000 m^2/kg
n_orbits=1000;
tspan=linspace(0,n_orbits*T,n_orbits*1000); 
perturbation_ID=[1,3];
odefun= @(t,y) ode_2bp_perturbed(t,y,perturbation_ID,parameters);
[t,y_perturbed] = ode113( odefun, tspan, y0, options);
%in the following 4 lines we are trasforming the integrated cartesian
%coodinates in to keplerian one
A=zeros(length(t),6);
for ii =1:length(t)
    [a,e,i,OM,om,th]=car2kep(y_perturbed(ii,1:3),y_perturbed(ii,4:6),mu_E);
    A(ii,:)=[a,e,i,OM,om,th];
end

figure()
Earth_3D(1)
hold on
patch([y_perturbed(:,1);nan]./R_e, [y_perturbed(:,2);nan]./R_e, [y_perturbed(:,3);nan]./R_e, [t;t(end)]./T, edgecolor='interp', facecolor='None')
H=colorbar;
xlabel('x [R_e]');
ylabel('y [R_e]');
zlabel('z [R_e]');
ylabel(H, 'orbit period [T]')


%REPEATING GROUND TRACK WITH DISTURBANCES 
% ground track with disturbances J2 & SRP uncorrected
[~,~] =groundtrack(t,y_perturbed,omega_e);

%same for repiting groundtrack
[t,y_disturbed_repeating]=ode113(odefun,tspan,y0_repeating_GT, options);

% ground track with disturbances J2 & SRP 
[lambda,phi] =groundtrack(t,y_disturbed_repeating,omega_e);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ASSIGNMENT 3
% INITIAL ORBIT PERTURBED IN KEPLERIAN COORDINATES

perturbation_ID = [2,4]; % disturbless, to see more details type 'help compute_disturb'.
odefun= @(t,y) gauss_ode_perturbed(t, y, perturbation_ID, parameters);
[~,x_kep__with_J2_SRP]=ode113(odefun,tspan,x0,options);

perturbation_ID = 0; % disturbless, to see more details type 'help compute_disturb'.
odefun= @(t,y) gauss_ode_perturbed(t, y, perturbation_ID, parameters);
[~,x_kep]=ode113(odefun,tspan,x0,options);

perturbation_ID = 2; % disturbless, to see more details type 'help compute_disturb'.
odefun= @(t,y) gauss_ode_perturbed(t, y, perturbation_ID, parameters);
[~,x_kep__with_J2]=ode113(odefun,tspan,x0,options);

perturbation_ID = 4; % disturbless, to see more details type 'help compute_disturb'.
odefun= @(t,y) gauss_ode_perturbed(t, y, perturbation_ID, parameters);
[t,x_kep_with_SRP]=ode113(odefun,tspan,x0,options);
% %  % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% PLOT 2: Keplerian element

% String to put in the legend
str = {'Earth';'Nominal orbit';'J2 perturbed orbit';
    'SRP perturbed orbit'; 'J2 and SRP perturbed orbit'};
% First part --> a e i;
fh = figure; fsize = 10;
title('x_kep -->sma e i')
subplot(2,3,1);
plot(t/T,x_kep(:,1)*1e-3,'r','Linewidth',2) %initial sma
hold on
plot(t/T,x_kep__with_J2(:,1)*1e-3,'g:','Linewidth',1.3) %perturbed sma in keplerian coordinates
plot(t/T,x_kep_with_SRP(:,1)*1e-3,'b:','Linewidth',1.3) %perturbed sma in keplerian coordinates
plot(t/T,x_kep__with_J2_SRP(:,1)*1e-3,'m:','Linewidth',1.3) %perturbed sma in keplerian coordinates
xlabel('Time [Orbits]','Interpreter','latex'),ylabel('SMA [$10^3km$]','Interpreter','latex')
ax = gca; ax.FontSize = fsize;

subplot(2,3,2);%
plot(t/T,x_kep(:,2),'r','Linewidth',2) %initial e
hold on
plot(t/T,x_kep__with_J2(:,2),'g:','Linewidth',1.3)
plot(t/T,x_kep__with_J2_SRP(:,2),'m:','Linewidth',1.3)
plot(t/T,x_kep_with_SRP(:,2),'b:','Linewidth',2)
xlabel('Time [Orbits]','Interpreter','latex'),ylabel('e [-]','Interpreter','latex')
title('Semi major axis, exentricity and inclination of the orbit perturbed, with SRP & J2 acceleration, and not  in Keplerian coordinates')
ax = gca; ax.FontSize = fsize;

subplot(2,3,3);
plot(t/T,wrapTo360(x_kep(:,3)*180/pi),'r','Linewidth',2) %initial e
hold on
plot(t/T,wrapTo360(x_kep__with_J2(:,3)*180/pi),'g:','Linewidth',1.3)
plot(t/T,wrapTo360(x_kep_with_SRP(:,3)*180/pi),'b:','Linewidth',1.3)
plot(t/T,wrapTo360(x_kep__with_J2_SRP(:,3)*180/pi),'m:','Linewidth',1.3)
xlabel('Time [Orbits]','Interpreter','latex'),ylabel('i [deg]','Interpreter','latex')
legend(str(2:5),'Interpreter','latex')
ax = gca; ax.FontSize = fsize;
fh.WindowState = 'maximized';


% Second part --> Om om theta
subplot(2,3,4);
plot(t/T,wrapTo360(x_kep(:,4)*180/pi),'r','Linewidth',2) %initial OM
hold on
plot(t/T,wrapTo360(x_kep__with_J2(:,4)*180/pi),'g:','Linewidth',1.3)
plot(t/T,wrapTo360(x_kep_with_SRP(:,4)*180/pi),'b:','Linewidth',1.3)
plot(t/T,wrapTo360(x_kep__with_J2_SRP(:,4)*180/pi),'m:','Linewidth',1.3)
xlabel('Time [Orbits]','Interpreter','latex'),ylabel('$\Omega$ [deg]','Interpreter','latex')
ax = gca; ax.FontSize = fsize;

subplot(2,3,5);
plot(t/T,wrapTo360(x_kep(:,5)*180/pi),'r','Linewidth',2) %initial om
hold on
plot(t/T,wrapTo360(x_kep__with_J2(:,5)*180/pi),'g:','Linewidth',1.3)
plot(t/T,wrapTo360(x_kep_with_SRP(:,5)*180/pi),'b:','Linewidth',1.3)
plot(t/T,wrapTo360(x_kep__with_J2_SRP(:,5)*180/pi),'m:','Linewidth',1.3)
xlabel('Time [Orbits]','Interpreter','latex'),ylabel('$\omega$ [deg]','Interpreter','latex')
title('Omega, omega and theta of the orbit perturbed, with SRP & J2 acceleration, and not  in Keplerian coordinates');
ax = gca; ax.FontSize = fsize;

subplot(2,3,6);
plot(t/T,wrapTo360(x_kep(:,6)*180/pi),'r','Linewidth',2) 
hold on
plot(t/T,wrapTo360(x_kep__with_J2_SRP(:,6)*180/pi),'m:','Linewidth',1.3)
xlabel('Time [Orbits]','Interpreter','latex'),ylabel('$\theta$ [deg]','Interpreter','latex')
ax = gca; ax.FontSize = fsize;
fh.WindowState = 'maximized';


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% % First part --> a e i;
% fh = figure; fsize = 10;
% subplot(1,3,1);
% plot(t/T,A(:,1)*1e-3,'r','Linewidth',2) %initial sma
% hold on
% plot(t/T,A(:,1)*1e-3,'m:','Linewidth',1.3) %perturbed sma in keplerian coordinates
% title('A -->sma e i');
% xlabel('Time [Orbits]','Interpreter','latex'),ylabel('SMA [$10^3km$]','Interpreter','latex')
% ax = gca; ax.FontSize = fsize;
% %il semimajoraxis diminuisce e poi riamunta perchè l'orbita ruora e poi aumenta il suo sma senza però
% 
% subplot(1,3,2);%variare l'energia dell'orbita stessa infatti srp non diminuisce l'energia
% plot(t/T,A(:,2),'r','Linewidth',2) %initial e
% hold on
% plot(t/T,A(:,2),'m:','Linewidth',1.3)
% xlabel('Time [Orbits]','Interpreter','latex'),ylabel('e [-]','Interpreter','latex')
% ax = gca; ax.FontSize = fsize;
% 
% subplot(1,3,3);
% plot(t/T,wrapTo360(A(:,3)*180/pi),'r','Linewidth',2) %initial e
% hold on
% plot(t/T,wrapTo360(A(:,3)*180/pi),'m:','Linewidth',1.3)
% xlabel('Time [Orbits]','Interpreter','latex'),ylabel('i [deg]','Interpreter','latex')
% legend(str([2:5]),'Interpreter','latex')
% ax = gca; ax.FontSize = fsize;
% fh.WindowState = 'maximized';
% %i cambi di inclinazione sono ascrivibili a j2
% 
% % Second part --> Om om theta;
% fh = figure; fsize = 10;
% subplot(1,3,1);%Omega scende con j2 quando l'inclinazione è minore di 90
% title('A -->Om om theta');
% plot(t/T,wrapTo360(A(:,4)*180/pi),'r','Linewidth',2) %initial OM
% hold on 
% plot(t/T,wrapTo360(A(:,4)*180/pi),'m:','Linewidth',1.3)
% xlabel('Time [Orbits]','Interpreter','latex'),ylabel('$\Omega$ [deg]','Interpreter','latex')
% ax = gca; ax.FontSize = fsize;
% 
% subplot(1,3,2);
% plot(t/T,wrapTo360(A(:,5)*180/pi),'r','Linewidth',2) %initial om
% hold on
% plot(t/T,wrapTo360(A(:,5)*180/pi),'m:','Linewidth',1.3)
% xlabel('Time [Orbits]','Interpreter','latex'),ylabel('$\omega$ [deg]','Interpreter','latex')
% legend(str([2:5]),'Interpreter','latex')
% ax = gca; ax.FontSize = fsize;
% 
% subplot(1,3,3);
% plot(t/T,A(:,6)*180/pi,'r','Linewidth',2) %initial the
% hold on
% plot(t/T,A(:,6)*180/pi,'m:','Linewidth',1.3)
% xlabel('Time [Orbits]','Interpreter','latex'),ylabel('$\theta$ [deg]','Interpreter','latex')
% ax = gca; ax.FontSize = fsize;
% fh.WindowState = 'maximized';
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%4.4 ERROR PLOT X-A WITH SEMILOGX
for ii =3:6
    x_kep__with_J2_SRP (:,ii)=wrapTo2Pi(x_kep__with_J2_SRP(:,ii));
end
    err=abs(A(:,1)-x_kep__with_J2_SRP(:,1))./sma;
    err2=abs(A(:,2)-x_kep__with_J2_SRP(:,2))./e;
    err3=abs(A(:,3)-x_kep__with_J2_SRP(:,3))./360;
    err4=abs(A(:,4)-x_kep__with_J2_SRP(:,4))./360;
    err5=abs(A(:,5)-x_kep__with_J2_SRP(:,5))./360;
    err6=abs(A(:,6)-x_kep__with_J2_SRP(:,6))./360;
%%

figure()
subplot(2,3,1);
semilogy(t/T,err);
title('a');
xlabel('time [T]');
ylabel('|a_{gauss}-a_{cart}|/ a_0 [km]');
subplot(2,3,2);
semilogy(t/T,err2);
title('e');
xlabel('time [T]');
ylabel('|e_{gauss}-e_{cart}|');
subplot(2,3,3);
semilogy(t/T,err3);
title('i');
xlabel('time [T]');
ylabel('|i_{gauss}-i_{cart}|/ 2\pi [deg]');
subplot(2,3,4);
semilogy(t/T,err4);
title('\Omega');
xlabel('time [T]');
ylabel('|\Omega_{gauss}-\Omega_{cart}|/ 2\pi  [deg]');
subplot(2,3,5);
semilogy(t/T,err5);
title('\Omega');
xlabel('time [T]');
ylabel('|\omega_{gauss}-\omega_{cart}|/ 2\pi [deg]');
subplot(2,3,6);
semilogy(t/T,err6);
title('\theta');
xlabel('time [T]');
ylabel('|\theta_{gauss}-\theta_{cart}|/\theta_{gauss}  [deg]');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Filtering
x_kep_with_J2_SRP(:, 1:2) = (x_kep__with_J2_SRP(:,1:2));
x_kep_with_J2_SRP(:, 3:6) = [rad2deg(x_kep__with_J2_SRP(:,3:6))];
A(:, 3:6) = [rad2deg(A(:,3:6))];
kep_f = movmean(A, 1000);
kep_f_2 = movmean(A, 80);

fh = figure();
subplot(2,3,1);
fsize = 10;
hold on
plot(t/T,x_kep_with_J2_SRP(:,1)*1e-3,'m:','Linewidth',1.3,'DisplayName','unfiltered data') %perturbed sma in keplerian coordinates
plot(t/T, kep_f_2(:, 1)*1e-3, 'g.', 'DisplayName', 'lowpass');
plot(t/T, kep_f(:, 1)*1e-3, 'b.', 'DisplayName', 'Secular');
title('a');
xlabel('time [T]');
ylabel('a [km]');ax = gca; ax.FontSize = fsize;


axes('position', [0.2 0.66 0.1 0.1])
hold on
plot(t/T,x_kep_with_J2_SRP(:,1)*1e-3,'m:','Linewidth',1.3) %perturbed sma in keplerian coordinates
plot(t/T, kep_f_2(:, 1)*1e-3, 'g.', 'DisplayName', 'lowpass');
plot(t/T, kep_f(:, 1)*1e-3, 'b.', 'DisplayName', 'Secular');
axis([t(100)/T, t(10000)/T, 18.8, 18.9])
fh.WindowState = 'maximized';



subplot(2,3,2);
hold on;
grid on;
plot(t/T, x_kep_with_J2_SRP(:, 2),'m:','Linewidth',1.3,'DisplayName','unfiltered data') %perturbed e in keplerian coordinates
plot(t/T, kep_f_2(:, 2), 'g.', 'DisplayName', 'lowpass');
plot(t/T, kep_f(:, 2), 'b.', 'DisplayName', 'Secular');
title('e');
xlabel('time [T]');
ylabel('e [-]');


axes('position', [0.45 0.66 0.1 0.1])
hold on
plot(t/T, x_kep_with_J2_SRP(:, 2),'m:','Linewidth',1.3);
plot(t/T, kep_f_2(:, 2),'g.','Linewidth',1.3);
plot(t/T, kep_f(:, 2), 'b.','Linewidth',1.3);
axis([t(100)/T, t(10000)/T, 0.578, 0.58])
fh.WindowState = 'maximized';



subplot(2,3,3);
hold on;
grid on;
plot(t/T, x_kep_with_J2_SRP(:, 3),'m:','Linewidth',1.3,'DisplayName','unfiltered data') %perturbed i in keplerian coordinates
plot(t/T, kep_f_2(:, 3), 'g.', 'DisplayName', 'lowpass');
plot(t/T, kep_f(:, 3), 'b.', 'DisplayName', 'Secular');
title('i');
xlabel('time [T]');
ylabel('i [deg]');
legend;

axes('position', [0.75 0.66 0.1 0.1])
hold on
plot(t/T, x_kep_with_J2_SRP(:, 3),'m:','Linewidth',1.3,'DisplayName','unfiltered data') %perturbed i in keplerian coordinates
plot(t/T, kep_f_2(:, 3), 'g.','Linewidth',1.3);
plot(t/T, kep_f(:, 3), 'b.','Linewidth',1.3);
axis([t(100)/T, t(10000)/T, 64.25, 64.30])
fh.WindowState = 'maximized';


subplot(2,3,4);
hold on;
grid on;
plot(t/T, x_kep_with_J2_SRP(:, 4),'m:','Linewidth',1.3,'DisplayName','unfiltered data') %perturbed Om in keplerian coordinates
plot(t/T, kep_f_2(:,4), 'g.', 'DisplayName', 'lowpass');
plot(t/T, kep_f(:, 4), 'b.', 'DisplayName', 'Secular');
title('\Omega');
xlabel('time [T]');
ylabel('\Omega [deg]');


axes('position', [0.2 0.23 0.1 0.1])
hold on
plot(t/T, x_kep_with_J2_SRP(:,4),'m:','Linewidth',1.3,'DisplayName','unfiltered data') %perturbed i in keplerian coordinates
plot(t/T, kep_f_2(:, 4), 'g.','Linewidth',1.3);
plot(t/T, kep_f(:,4), 'b.','Linewidth',1.3);
axis([t(100)/T, t(10000)/T, 179.87, 180])
fh.WindowState = 'maximized';


subplot(2,3,5);
hold on;
grid on;
plot(t/T, x_kep_with_J2_SRP(:, 5),'m:','Linewidth',1.3,'DisplayName','unfiltered data') %perturbed Om in keplerian coordinates
plot(t/T, kep_f_2(:,5), 'g.', 'DisplayName', 'lowpass');
plot(t/T, kep_f(:, 5), 'b.', 'DisplayName', 'Secular');
title('\omega');
xlabel('time [T]');
ylabel('\omega [deg]');


axes('position', [0.45 0.23 0.1 0.1])
hold on
plot(t/T, x_kep_with_J2_SRP(:, 5),'m:','Linewidth',1.3,'DisplayName','unfiltered data') %perturbed i in keplerian coordinates
plot(t/T, kep_f_2(:, 5), 'g.','Linewidth',1.3);
plot(t/T, kep_f(:, 5), 'b.','Linewidth',1.3);
axis([t(100)/T, t(10000)/T, 179.87, 180])
fh.WindowState = 'maximized';

%% Check computational time

sim_size=10;
n_orbits=10;
tspan=linspace(0, n_orbits*T, n_orbits*1000); 

simulation_2bp=zeros(sim_size,1);
simulation_gauss=zeros(sim_size,1);
for jj=1:sim_size
    perturbation_ID = [1,3]; % disturbless, to see more details type 'help compute_disturb'.
    odefun= @(t,y) ode_2bp_perturbed(t, y, perturbation_ID, parameters);
    tic()
    [t,y_repeating_GT]=ode113(odefun,tspan,y0,options);
    simulation_2bp(jj)=toc();

    perturbation_ID = [2,4]; % disturbless, to see more details type 'help compute_disturb'.
    odefun= @(t,y) gauss_ode_perturbed(t, y, perturbation_ID, parameters);
    tic()
    [~,x_kep__with_J2_SRP]=ode113(odefun,tspan,x0,options);
    simulation_gauss(jj)=toc();
end

figure;
plot(simulation_2bp)
hold on
plot(simulation_gauss)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Horizons results from COSMOS-2546
%Start date: 9 May 2020 00:00
start_date_mjd2000 = 58978 - 51544.5;
end_date_mjd2000 = 59726 - 51544.5;
%End date: 27 May 2022 00:00
%   Data point every 8 hours
A=importdata('ephemerides.txt');
e_sc=A.data(:,1);%eccentricity
%r_p_sc=A.data(:,2); %periaxis [km]
i_sc=A.data(:,3).*pi/180;%inclination [rad]
Omega_sc=A.data(:,4).*pi/180;%ascending node [rad]
% Omega_sc=pi*ones(length(i),1); [rad]
omega_sc=A.data(:,5).*pi/180;%argument of pericenter [rad]
% omega_sc=pi*ones(length(i),1); [rad]
%T_p_sc=A.data(:,6);%periaxis time [s]
%n_sc=A.data(:,7);%mid motion
%M_sc=A.data(:,8);%M angle
theta_sc=A.data(:,9).*pi/180;%true anomaly [rad]
sma_sc=A.data(:,10);%semimajor axis [km]
%r_a_sc=A.data(:,11);%apoaxis [km]
%T_sc=A.data(:,12);%sidereal period [s]
T=2*pi*sqrt(sma_sc(1000)^3/mu_E); %sec %period of the sc around earth
kep0= [sma_sc(1),e_sc(1), i_sc(1), Omega_sc(1),omega_sc(1), theta_sc(1)];

T_max=length(sma_sc)*8*3600;
tspan=linspace(0,T_max, 1000000);

t_sc= linspace(0,2*365.25*24*3600,length(sma_sc));

%%
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14);
perturbation_ID = [2,4]; % disturbless, to see more details type 'help compute_disturb'.
odefun= @(t,y) gauss_ode_perturbed(t, y, perturbation_ID, parameters);
[~,x_kep]=ode113(odefun,tspan,kep0,options);

%% PLOT 

% String to put in the legend
str = {'Earth';'Nominal orbit';'J2 perturbed orbit';
    'SRP perturbed orbit'; 'J2 and SRP perturbed orbit'};
% First part --> a e i;
fh = figure; fsize = 10;
title('x_kep -->sma e i')
subplot(1,3,1);
plot(tspan/T,x_kep(:,1)*1e-3,'b','Linewidth',1.3) 
hold on
plot(t_sc/T,sma_sc(:,1)*1e-3,'r:','Linewidth',1.3) 
xlabel('Time [Orbits]','Interpreter','latex'),ylabel('SMA [$10^3km$]','Interpreter','latex')
ax = gca; ax.FontSize = fsize;

subplot(1,3,2);
plot(tspan/T,x_kep(:,2),'b','Linewidth',1.3)
hold on
plot(t_sc/T,e_sc(:),'r:','Linewidth',1.3)
xlabel('Time [Orbits]','Interpreter','latex'),ylabel('e [-]','Interpreter','latex')
title('Semi major axis, exentricity and inclination of the orbit perturbed, with SRP & J2 acceleration, and not  in Keplerian coordinates')
ax = gca; ax.FontSize = fsize;

subplot(1,3,3);
plot(tspan/T,wrapTo360(x_kep(:,3)*180/pi),'b','Linewidth',1.3) 
hold on
plot(t_sc/T,wrapTo360(i_sc(:)*180/pi),'r:','Linewidth',1.3)
xlabel('Time [Orbits]','Interpreter','latex'),ylabel('i [deg]','Interpreter','latex')
legend(str(2:5),'Interpreter','latex')
ax = gca; ax.FontSize = fsize;
fh.WindowState = 'maximized';


% Second part --> Om om theta
fh = figure; fsize = 10;

subplot(1,2,1);
plot(tspan/T,wrapTo360(x_kep(:,4)*180/pi),'b','Linewidth',1.3)
hold on
plot(t_sc/T,wrapTo360(Omega_sc(:)*180/pi),'r:','Linewidth',1.3)
xlabel('Time [Orbits]','Interpreter','latex'),ylabel('$\Omega$ [deg]','Interpreter','latex')
ax = gca; ax.FontSize = fsize;

omega_f= movmean(x_kep(:,5),1000000);

subplot(1,2,2)
plot(tspan/T,wrapTo360(x_kep(:,5)*180/pi),'b','Linewidth',1.3) %simulation result
hold on
plot(t_sc/T,wrapTo360(omega_sc(:)*180/pi),'r:','Linewidth',1.3)
plot(tspan/T,wrapTo360(omega_f(:)*180/pi),'g:','Linewidth',1.3)
xlabel('Time [Orbits]','Interpreter','latex'),ylabel('$\omega$ [deg]','Interpreter','latex')
title('Omega, omega and theta of the orbit perturbed, with SRP & J2 acceleration, and not  in Keplerian coordinates');
legend(str(2:5),'Interpreter','latex')
ax = gca; ax.FontSize = fsize;
