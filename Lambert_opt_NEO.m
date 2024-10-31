function  [delta_v_3,v_initial_Lambert_3_x,v_initial_Lambert_3_y,v_initial_Lambert_3_z] = Lambert_opt_NEO(mjd_departure_in,mjd_departure_fin,mjd_arrival_in,mjd_arrival_fin,unit,Planet_1,NEO,orbitType,delta_v_1_MAX)

% MOST EFFICIENT ( MIN delta_v) LAMBERT TRANSFER BETWEEN A STARTING PLANET
% (Planet1) AND AN ARRIVAL NEO (NEO) GIVEN A DEPARTURE WINDOW AND AN
% ARRIVAL WINDOW COMPUTING ONLY delta_v_3 TO GO FROM THE LAMBERT LEG TO THE
% NEO's ORBIT

%NOTE. this function uses only a zero revolution transfer 
%NOTE: this function gives the CONSTRAIN on the Maximum excess velocity from launcher (constrain on delta_v1)


% PLANET'S LEGEND:
%                   1:   Mercury
%                   2:   Venus
%                   3:   Earth
%                   4:   Mars
%                   5:   Jupiter
%                   6:   Saturn
%                   7:   Uranus
%                   8:   Neptune
%                   9:   Pluto
%                   10:  Sun


% INPUT 
 % departure window in mjd: mjd_departure_in,mjd_departure_fin
 % arrival window in mjd: mjd_arrival_in,mjd_arrival_fin
 % time step between a point and the next in departure and arrival
 % windows , expressed in days : unit

 % the number related to starting planet : Planet1
 % the number related to arrival NEO : NEO
 % the options used in the Lambert solver:orbitType    Logical variable defining whether transfer is
 %                       0: direct transfer from R1 to R2 (counterclockwise)
 %                       1: retrograde transfer from R1 to R2 (clockwise)
 % the constrain on delta_v_3: delta_v_1_MAX


% OUTPUT
 % the delta_v_3 needed for the transfer for each departure point and each arrival point: delta_v_3
 % the components of the velocity at the start of the second Lambert arc: v_initial_Lambert_3_x
 %                                                                        v_initial_Lambert_3_y
 %                                                                        v_initial_Lambert_3_z

 
% CONTRIBUTORS
 % Francesco Paolo Vacca



% MODIFIED JULIAN DAYS
% DEPARTURE WINDOW

mjd_departure=[mjd_departure_in:unit:mjd_departure_fin];
% ARRIVAL WINDOW

mjd_arrival=[mjd_arrival_in:unit:mjd_arrival_fin];

% POSITION AND VELOCITY OF THE SPACECRAFT ON THE INITIAL ORBIT IN THE MOMENT OF DEPARTURE 

mjd2000_departure=mjd_departure - 51544.5; % mjd2000 day for departure
ibody_1_vect=Planet_1*ones(1,length(mjd2000_departure));
for i=1:length(mjd2000_departure)
[kep_1(i,:),ksun] = uplanet(mjd2000_departure(i), ibody_1_vect(i));  %orbital parameters wrt the starting point for orbit propagation
a_1(i)=kep_1(i,1);
e_1(i)=kep_1(i,2);
i_1(i)=kep_1(i,3);
OM_1(i)=kep_1(i,4);
om_1(i)=kep_1(i,5);
theta_0_1(i)=kep_1(i,6);

% initial cartesian parameters
[rr_0_1(i,:),vv_0_1(i,:)] = par2car(a_1(i),e_1(i),i_1(i),OM_1(i),om_1(i), theta_0_1(i),ksun);
end

% POSITION AND VELOCITY OF THE SPACECRAFT ON THE FINAL ORBIT IN THE MOMENT OF ARRIVAL

mjd2000_arrival=mjd_arrival - 51544.5; % mjd2000 day for arrival
ibody_2_vect=NEO*ones(1,length(mjd2000_arrival));
for j=1:length(mjd2000_arrival)
[kep_NEO(j,:),mass_NEO(j),M(j)] = ephNEO(mjd2000_arrival(j),NEO); %orbital parameters wrt the arriving point for orbit propagation
a_2(j)=kep_NEO(j,1);
e_2(j)=kep_NEO(j,2);
i_2(j)=kep_NEO(j,3);
OM_2(j)=kep_NEO(j,4);
om_2(j)=kep_NEO(j,5);
theta_0_2(j)=kep_NEO(j,6);

% initial cartesian parameters
[rr_0_2(j,:),vv_0_2(j,:)] = par2car(a_2(j),e_2(j),i_2(j),OM_2(j),om_2(j), theta_0_2(j),ksun);
end

% TRANSFER ORBIT
% solve Lambert problem for the 2 given orbits, wrt the departure point
% and the arrival point. Compute delta_v needed to go on the transfer orbit
% from initial orbit on the departure point and compute delta_v needed to go
% from fransfer obit to final orbit on the arrival point

Nrev=0; %number of revolutions
optionsLMR=1;%warnings are displayed only when the algorithm does not converge


for i=1:length(mjd2000_departure)
    for j=1:length(mjd2000_arrival)
        delta_t_matrix(i,j)=(mjd2000_arrival(j) - mjd2000_departure(i))*24*60*60 ;% TIME OF FLIGTH [s]
%following the column departure day changes, arrival day is remains the same
%following the row arrival day changes, departure day is remains the same

%starting from each departing point (i) and arriving in each arrival
%point (j)

       [a_T_matr(i,j),p_T_matr(i,j),e_T_matr(i,j),ERROR_matr(i,j),v_i_matr(i,j,:),v_f_matr(i,j,:),TPAR_matr(i,j),THETA_matr(i,j)] = lambertMR(rr_0_1(i,:),rr_0_2(j,:),delta_t_matrix(i,j),ksun,orbitType,Nrev,optionsLMR);

 v_i_matr(abs(v_i_matr)>100) = NaN;

       v_initial_Lambert_3_x=v_i_matr(:,:,1); 
       v_initial_Lambert_3_y=v_i_matr(:,:,2); 
       v_initial_Lambert_3_z=v_i_matr(:,:,3); 

       v_f_x=v_f_matr(:,:,1); 
       v_f_y=v_f_matr(:,:,2); 
       v_f_z=v_f_matr(:,:,3); 


       %components of the delta_v3 matrix, delta_v needed to go from
       %transfer orbit to second orbit changhing departure and arrival point
       delta_v_3_x(i,j)=vv_0_2(j,1) - v_f_x(i,j);
       delta_v_3_y(i,j)=vv_0_2(j,2) - v_f_y(i,j);
       delta_v_3_z(i,j)=vv_0_2(j,3) - v_f_z(i,j);

       
       delta_v_3(i,j)=sqrt(delta_v_3_x(i,j)^2 + delta_v_3_y(i,j)^2 + delta_v_3_z(i,j)^2);

       delta_v_tot(i,j)=delta_v_3(i,j);   %total delta_v needed for the mission, each column has a fixed arrival day and a changing departure day.
       %each row has a fixed departure day and a changing arrival day.

       if delta_v_3(i,j) > delta_v_1_MAX
           delta_v_3(i,j)=NaN;
       end

    end
end
 

% figure
% surf(datenum(mjd2000_departure),datenum(mjd2000_arrival),delta_v_tot')
% datetick('x','dd-mmm-yy'),datetick('y','dd-mmm-yy');
% colorbar
% xlabel 'Departure date'
% ylabel 'Arrival date'
% zlabel 'delta_v [km/s]'


% figure
% plot(datenum(mjd2000_arrival),delta_v_tot')
% datetick('x','dd-mmm-yy');
% 
% xlabel 'Arrival date'
% ylabel 'delta_v [Km/s]'



% figure
% contour(datenum(mjd2000_departure),datenum(mjd2000_arrival),delta_v_tot',ShowText="on")
% datetick('x','dd-mmm-yy'),datetick('y','dd-mmm-yy');
% colorbar 
%  xlabel 'Departure date'
% ylabel 'Arrival date'
% zlabel 'delta_v [km/s]'
% grid on