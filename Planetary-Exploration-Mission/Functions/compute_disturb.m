function [a_tot,a_vec]= compute_disturb(t, y, perturbation_ID, parameters)
% perturbation_ID is the vector containing the one or multiple disturb
% among:
% perturbation_ID = 0 or > 4 --> No disturb
% perturbation_ID = 1 --> J2 acceleration in carthesian coordinates
% perturbation_ID = 2 --> J2 acceleration in Keplerian coordinates
% perturbation_ID = 3 --> SRP acceleration in carthesian coordinates
% perturbation_ID = 4 --> SRP acceleration in Kepleriane coordinates

% AUTHOR: 
 % Lorenzo Dionigi

% SUPERVISOR:
 % Prof. Camilla Colombo


% Opening the inputs

R_e = parameters.R_e;
AU = parameters.AU;
mu = parameters.mu;
J2 = parameters.J2;
Cr = parameters.Cr;
A_M = parameters.A_M;
p_Solar = parameters.p_Solar;
% Main cycle
a_tot = zeros(3,1);
for ii = 1:length(perturbation_ID) % It runs all the disturbs
    switch perturbation_ID(ii)
        case 0 % No disturb is applied
            a = 0;

        case 1 % J2 disturb in cathesian coordinates
            r = y(1:3);
            rnorm = norm(r);
            coeff = 3/2*J2*mu*R_e^2/rnorm^4;
            a_x = coeff*(r(1)/rnorm*(5*r(3)^2/rnorm^2 -1));
            a_y = coeff*(r(2)/rnorm*(5*r(3)^2/rnorm^2 -1));
            a_z = coeff*(r(3)/rnorm*(5*r(3)^2/rnorm^2 -3));
            a = [a_x;a_y;a_z];

        case 2 % J2 disturb in Keplerian coordinates (RSW)
            a = y(1);
            e = y(2);
            i = y(3);
            % OM = y(4); Not needed
            om = y(5);
            th = y(6);

            p = a*(1-e^2);
            r = p/(1 + e*cos(th));

            coeff = -3/2*J2*mu*R_e^2/r^4;
            a_r = coeff*(1 - 3*sin(i)^2*sin(th + om)^2);
            a_s = coeff*sin(i)^2*sin(2*(th+om));
            a_w = coeff*sin(2*i)*sin(th+om);

            a = [a_r; a_s; a_w];

        case 3 % SRP disturb in carthesian
            r = y(1:3);
            Date0 =  date2mjd2000([2020,01,01,00,00,00]); % reference time to start simulation in mjd2000
            Date_actual = Date0 + t/3600/24; % In Days
            [kep_Earth, mu_Sun] = uplanet(Date_actual, 3);
            R_earth = kep2car(kep_Earth(1),kep_Earth(2),kep_Earth(3),kep_Earth(4),kep_Earth(5),kep_Earth(6),mu_Sun); % This is the Earth position wrt Sun
            R_SUN = - R_earth;  % Thus, this is the Sun position wrt Earth
            r_SC_SUN = r - R_SUN; % Vector of distance from Spacecraft to Sun
            rnorm_SC_SUN = norm(r_SC_SUN);
            coeff = p_Solar*AU^2/(rnorm_SC_SUN^2)*Cr*A_M; % m/s^2
            a = -coeff*1e-3*r_SC_SUN/rnorm_SC_SUN; % km/s^2

        case 4 % SRP disturb in Keplerian
            [r,v] = kep2car(y(1),y(2),y(3),y(4),y(5),y(6),mu);
            Date0 =  date2mjd2000([2020,01,01,00,00,00]); % reference time to start simulation in mjd2000
            Date_actual = Date0 + t/3600/24; % In Days
            [kep_Earth, mu_Sun] = uplanet (Date_actual, 3);
            R_earth = kep2car(kep_Earth(1),kep_Earth(2),kep_Earth(3),kep_Earth(4),kep_Earth(5),kep_Earth(6),mu_Sun); % This is the Earth position wrt Sun
            R_SUN = - R_earth;  % Thus, this is the Sun position wrt Earth
            r_SC_SUN = r - R_SUN; % Vector of distance from Spacecraft to Sun
            rnorm_SC_SUN = norm(r_SC_SUN);
            coeff = p_Solar*AU^2/(rnorm_SC_SUN^2)*Cr*A_M; % m/s^2
            a = -coeff*1e-3*r_SC_SUN/rnorm_SC_SUN; % km/s^2 in ECEI
            
            % Conversion of a_drag from ECEI -> RSW
            l1 = r/norm(r);
            h = cross(r,v);
            l3 = h/norm(h);
            l2 = cross(l3,l1);
            ALN = [l1, l2, l3]';

            a = ALN*a; % Rotation from EXEI to RSW frame

        otherwise
            a = 0;
    end % End of the 'switch'
    a_tot = a_tot + a;
end % End of the 'for'

end
