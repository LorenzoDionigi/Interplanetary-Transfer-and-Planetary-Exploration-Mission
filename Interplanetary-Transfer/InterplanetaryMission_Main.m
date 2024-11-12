%% EARTH-SATURN-NEO OPTIMAL TRANSFER
% FIND THE OPTIMAL TRANSFER IN TERMS OF DELTA_v WITH THE AIM OF GOING FROM EARTH TO SATURN FOLLOWING A FIRST LAMBERT LEG,
% PERFORM A FLYBY AROUND SATURN AND THEN GOING FROM SATURN TO NEO(88) FOLLOWING A SECOND LAMBERT LEG, WITH THE
% CONSTRAIN OF THE FIRST POSSIBLE DEPARTURE FROM EARTH ON 31/09/2031 AND LATEST ARRIVAL ON NEO ON 28/03/66.
% In order to find every possible solution to the problem, the aviable time window has been divided in smaller windows :
% the departure window from Earth has be decided to be equal to a synodic period between Earth and Saturn, expoiting the 
% geometrical pattern between the orbits wich repetes itself every sinodic period.
% Same has been applyed for the arrival windows on Saturn and for the arrival windows on NEO.


%% IN ORDER TO SEE EARTH-SATURN-NEO TRANSFER: CASE= 'FLYBY TRANSFER'
%% IN ORDER TO SEE DIRECT EARTH-NEO TRANSFER: CASE= 'DIRECT TRANSFER'
addpath('./Functions/');
CASE= 'FLYBY TRANSFER';


switch CASE 
    case 'FLYBY TRANSFER'


clc
clear
close all

% COMPUTATION OF THE EPHEMERIS OF EARTH,SATURN AND NEO FOR AN ARBITRARY DAY
% (31/09/2031) IN ORDER TO OBTAIN THE PERIODS OF THE ORBIT AND THE SINODIC
% PERIOD BETWEEN EARTH AND SATURN:

%EARTH ORBIT
mjd_departure_in=63139;
mjd2000_departure=mjd_departure_in - 51544.5;
ibody_E_vect=3;
[kep_E_in,ksun] = uplanet(mjd2000_departure, ibody_E_vect);
mu_Sun=astroConstants(4);
%NEO ORBIT
id=88;
[kep_NEO_in,mass_NEO,M] = ephNEO(mjd2000_departure,id);
%SATURN ORBIT
ibody_Sat_vect=6;
[kep_Sat_in,ksun] = uplanet(mjd2000_departure, ibody_Sat_vect);

Planet_1=3;
Planet_2=6;
NEO=88;
orbitType=0;
% Physical constrain on the max delta_v that the spacecraft is able to grant:
delta_v_1_MAX=25;
delta_v_3_MAX=25;
mjd_departure_in=63139;

n_ES=1;
T_E=2*pi*sqrt(kep_E_in(1).^3/mu_Sun);
T_Sat=2*pi*sqrt(kep_Sat_in(1).^3/mu_Sun);
T_syn_ES=((T_E*T_Sat)/abs(T_E-T_Sat))/(60*60*24);
T_NEO=2*pi*sqrt(kep_NEO_in(1).^3/mu_Sun);
T_syn_NEOS=((T_NEO*T_Sat)/abs(T_NEO-T_Sat))/(60*60*24);
mjd_departure_fin= mjd_departure_in + n_ES*T_syn_ES;

% time span from first possible departure day from Earth to last possible
% arrival day on NEO
MAX_time=75600 - 63139;
% number of departure windows in the aviable timespan
N_1=floor(MAX_time/(n_ES*T_syn_ES));
unit=30; %time discretization


%% STEP 0 : COUNTOUR PLOTS
% In order to visualize and study any possible pattern on the delta_v needed for the transfers Earth_Saturn and Saturn_NEO, has been performed
% a computation of delta_v_1 and delta_v_3 where delta_v_1 is the delta_v needed for the transfer from Earth's orbit to the first Lambert leg,
% delta_v_3 is the delta_v needed for the transfer from second Lambert leg to NEO's orbit. ATTENTION: this computation does not consider the flyby.
% In order to observe properly the results, delta_v_1 will be computed for several Earth-Saturn synodic periods (n_1) and delta_v_3 will be computed
% for several Saturn-NEO synodic periods (n_2). The time discretization used is unit_1 (days)

n_1=10;
n_2=3;
unit_1=10;


mjd_departure_fin_1=mjd_departure_in+n_1*T_syn_ES;
mjd_arrival_in_1=mjd_departure_fin_1;
mjd_arrival_fin_1=mjd_arrival_in_1 + (n_1+0.1)*T_syn_ES;

[delta_v_1_1,v_final_Lambert_1_x_1,v_final_Lambert_1_y_1,v_final_Lambert_1_z_1] = Lambert_opt(mjd_departure_in,mjd_departure_fin_1,mjd_arrival_in_1,mjd_arrival_fin_1,unit_1,Planet_1,Planet_2,orbitType,delta_v_1_MAX);

mjd2000_departure_1=[mjd_departure_in:unit_1:mjd_departure_fin_1] - 51544.5;
mjd2000_arrival_1=[mjd_arrival_in_1:unit_1:mjd_arrival_fin_1] - 51544.5;

figure
subplot(1,2,1)
contour(datenum(mjd2000_departure_1),datenum(mjd2000_arrival_1),delta_v_1_1',ShowText="off")
datetick('x','dd-mmm-yy'),datetick('y','dd-mmm-yy');
colorbar 
 xlabel 'Departure date from Earth'
ylabel 'Arrival date on Saturn'
zlabel 'delta_v_1 [km/s]'
grid on
title '\DeltaVdep'

mjd_departure_in_2=mjd_arrival_in_1;
mjd_departure_fin_2=mjd_departure_in_2+n_2*T_syn_NEOS;
mjd_arrival_in_2=mjd_departure_fin_2;
mjd_arrival_fin_2=mjd_arrival_in_2 + (n_2+0.1)*T_syn_NEOS;

[delta_v_3_2,v_final_Lambert_1_x_2,v_final_Lambert_1_y_2,v_final_Lambert_1_z_2] = Lambert_opt_NEO(mjd_departure_in_2,mjd_departure_fin_2,mjd_arrival_in_2,mjd_arrival_fin_2,unit_1,Planet_2,NEO,orbitType,delta_v_1_MAX);

mjd2000_departure_2=[mjd_departure_in_2:unit_1:mjd_departure_fin_2] - 51544.5;
mjd2000_arrival_2=[mjd_arrival_in_2:unit_1:mjd_arrival_fin_2] - 51544.5;

subplot(1,2,2)
contour(datenum(mjd2000_departure_2),datenum(mjd2000_arrival_2),delta_v_3_2',ShowText="off")
datetick('x','dd-mmm-yy'),datetick('y','dd-mmm-yy');
colorbar 
 xlabel 'Departure date from Saturn'
ylabel 'Arrival date on NEO'
zlabel 'delta_v_3 [km/s]'
grid on
title '\DeltaVarr'

%% STEP 1 : FIRST APPROXIMATE SOLUTION
% Find all the possible solutions for the problem (in terms of delta_v_1 needed to move from Earth's orbit to the first Lambert Leg, delta_v_flyby needed
% to move from the first to the second Lambert leg while exploiting a flyby of Saturn, and delta_v_3 needed to move from the second Lambert leg to NEO's orbit)
% for all the aviable departure days from Earth and all the aviable arrival days on Saturn, both divided in N_1 departure and arrival windows.
% Then, considering the flyby as an impulsive manouvre, and conseguently choosing the departure days from Saturn to NEO equal to the arrival days on Saturn from Earth, 
% the second Lambert Leg was computed for any possible arrival day on NEO (divided in N_1 arrival windows).
% In this first approximation, any point has been choose with a 30-days discretization in order to perform a first coarse grid search and find the
% best 30-days long regions where the solution is convenient.

% COMPUTATION OF delta_v_1 for any day of departure of any departure window (k) from Earth, and for any arrival day of any arrival window (l) on Saturn.

for k=1:N_1-1  %DEPARTURE WINDOWS FROM EARTH: built with initial and final departure days

    mjd_departure_in_vect(k)=mjd_departure_in + (k-1)*n_ES*T_syn_ES;
    mjd_departure_fin_vect(k)= mjd_departure_in_vect(k) + n_ES*T_syn_ES;

       for l=k:N_1-1  %ARRIVAL WINDOWS ON SATURN: built with initial and final arrival days

    mjd_arrival_in_vect(l)=mjd_departure_fin + (l-1)*n_ES*T_syn_ES;
    mjd_arrival_fin_vect(l)= mjd_arrival_in_vect(l) + n_ES*T_syn_ES;

       %ARRIVING ON SATURN IN EVERY ARRIVAL WINDOW DEPARTING FROM EARTH AT EVERY DEPARTURE WINDOW
    [delta_v_1(:,:,k,l),v_final_Lambert_1_x(:,:,k,l),v_final_Lambert_1_y(:,:,k,l),v_final_Lambert_1_z(:,:,k,l)] = Lambert_opt(mjd_departure_in_vect(k),mjd_departure_fin_vect(k),mjd_arrival_in_vect(l),mjd_arrival_fin_vect(l),unit,Planet_1,Planet_2,orbitType,delta_v_1_MAX);

   end
end

% COMPUTATION OF delta_v_3 for any day of departure of any departure window (h) from Saturn (equal to the arrival days on Saturn),
% and for any arrival day of any arrival window (j) on NEO.

for h=1: N_1-1 %DEPARTURE WINDOWS FROM SATURN
   for j = h: N_1-1 %ARRIVAL WINDOWS ON NEO
       for i=1:ceil(n_ES*T_syn_ES/unit) %departing day from Saturn (number of possible departure days given by the width of the departing window divided by unit of discretization)
           mjd_arrival_SAT_matrix(h,:)=mjd_arrival_in_vect(h):unit:mjd_arrival_in_vect(h) + n_ES*T_syn_ES; %every row rapresents the arrival day on Saturn for every arrival window (h)
           mjd_arrival_NEO_matrix(j,:)= (mjd_arrival_in_vect(j) + n_ES*T_syn_ES):unit:(mjd_arrival_in_vect(j) + 2*n_ES*T_syn_ES); %every row rapresents the arrival day on NEO for every arrival window (j)
 
           [delta_v_3(i,:,h,j),v_initial_Lambert_3_x(i,:,h,j),v_initial_Lambert_3_y(i,:,h,j),v_initial_Lambert_3_z(i,:,h,j)] = Lambert_opt_NEO(mjd_arrival_SAT_matrix(h,i),mjd_arrival_SAT_matrix(h,i),mjd_arrival_NEO_matrix(j,1),mjd_arrival_NEO_matrix(j,end),unit,Planet_2,NEO,orbitType,delta_v_1_MAX);
       end
   end
end


% COARSE GRID SEARCH: search of inaccurate 30-days-long discretized solutions, found using a tuning method wich exludes all the solutions with delta_v_1+delta_v_3 > delta_v_MAX. 
% delta_v_MAX has been found by minimazing the total cost (delta_v_tot=delta_v_1+delta_v_p+delta_v_3) for inaccurate solutions.
% Following this framework has been found out a solution whit delta_v_tot=16.5 [Km/s]. This value has been used as a preliminar constrain since, 
% the existence of a solution with delta_v_tot < delta_v_MAX implies that all the solutions with delta_1 + delta_v3 > delta_v_MAX are
% less convenient and then they can be discarded a priori.


delta_v_MAX=16.5;
z=1;
for I=1:size(delta_v_1,1) % departing day from Earth
    for J=1:size(delta_v_1,2) % arrival day on Saturn = departing day from saturn
        for K=1:size(delta_v_1,3) % departing window from Earth
            for L=K:size(delta_v_1,4) % arrival window on Saturn = departing window from saturn
                for P=1:size(delta_v_3,2) % arrival day on NEO
                    for Q=L:size(delta_v_3,3) % arrival window on NEO


                         delta_v_TOT_preliminar(I,J,K,L,P,Q)=delta_v_1(I,J,K,L) + delta_v_3(J,P,L,Q);
                        
                      
% collecting all the inaccurate solutions that respect the preliminary constrain:
                       if  delta_v_TOT_preliminar(I,J,K,L,P,Q)< delta_v_MAX
                           
                          delta_v_TOT_good_preliminar(z)= delta_v_TOT_preliminar(I,J,K,L,P,Q);
                          I_vect(z)=I;
                          J_vect(z)=J;
                          K_vect(z)=K;
                          L_vect(z)=L;
                          P_vect(z)=P;
                          Q_vect(z)=Q;
                         
                          
                          delta_v1_vect(z) = delta_v_1(I,J,K,L);
                        
                          delta_v_3_vect(z)=delta_v_3(J,P,L,Q);

                          z=z+1;
                       end

                    end
                end
            end
        end
    end
end

%%  STEP 2: ACCURATE SOLUTION
% Given the preliminary low-accurate solutions that respect the first constrain, a more accurate search of the solution has been accomplished:
% As a first step, the mjd relative to the first of the 30-days windows of the low-accurate solutions have been computed, and from them have been
% built 30-days-long windows with a time discretization of 1 day for each of the low-accurate solution.

for numero=1:length(delta_v_TOT_good_preliminar)

% initial days for low-accurate solutions
mjd_departure_MIN_EARTH=mjd_departure_in + (K_vect(numero)-1).*n_ES.*T_syn_ES + (I_vect(numero) -1).*unit;
mjd_arrival_MIN_SAT=(mjd_departure_in + n_ES.*T_syn_ES) + (L_vect(numero)-1).*n_ES.*T_syn_ES + (J_vect(numero) -1).*unit;
mjd_arrival_MIN_NEO=(mjd_departure_in + 2.*n_ES.*T_syn_ES) + (Q_vect(numero)-1).*n_ES.*T_syn_ES + (P_vect(numero) -1).*unit;


 unit_good=1;
 
% 30-days windows for each low-accurate solution
mjd_departure_EARTH_WINDOW=mjd_departure_MIN_EARTH:unit_good:(mjd_departure_MIN_EARTH+unit);
mjd_arrival_SAT_WINDOW=mjd_arrival_MIN_SAT:unit_good:(mjd_arrival_MIN_SAT+unit);
mjd_arrival_NEO_WINDOW=mjd_arrival_MIN_NEO:unit_good:(mjd_arrival_MIN_NEO+unit);


 [delta_v_1_good,v_final_Lambert_1_x_good,v_final_Lambert_1_y_good,v_final_Lambert_1_z_good] = Lambert_opt(mjd_departure_MIN_EARTH,mjd_departure_EARTH_WINDOW(end),mjd_arrival_MIN_SAT,mjd_arrival_SAT_WINDOW(end),unit_good,Planet_1,Planet_2,orbitType,delta_v_1_MAX);


for H=1:length(mjd_arrival_SAT_WINDOW) %departing day from Saturn in the good window
[delta_v_3_good(H,:),v_initial_Lambert_3_x_good(H,:),v_initial_Lambert_3_y_good(H,:),v_initial_Lambert_3_z_good(H,:)] = Lambert_opt_NEO(mjd_arrival_SAT_WINDOW(H),mjd_arrival_SAT_WINDOW(H),mjd_arrival_NEO_WINDOW(1),mjd_arrival_NEO_WINDOW(end),unit_good,Planet_2,NEO,orbitType,delta_v_3_MAX);
end


% given the windows for departure from Earth, arrival on Saturn and
% arrival on NEO, delta_v_tot depends on the day of departure from
% Earth (X), the day of departure/arrival from/to Saturn (Y) and from the
% day of arival on NEO (Z).


% REFINED GRID SEARCH: search of the accurate 1-day-long discretizated solutions collecting and then minimizing the solution found 
% in the previous regions in order to find the minimum delta_v_preliminary and the corrisponding days for each region.

w=1;
for X=1:size(delta_v_1_good,1) % departing day from Earth
    for Y=1:size(delta_v_1_good,2) % arrival day on Saturn = departing day from Saturn
        for Z=1:size(delta_v_3_good,2) % arrival day on NEO

             delta_v_preliminary_good(X,Y,Z)=delta_v_1_good(X,Y) + delta_v_3_good(Y,Z);

  %collect all the accurate solutions inside each window formed by
  %the inaccurate solutions
                        if  delta_v_preliminary_good(X,Y,Z)< delta_v_MAX 
                           
                          delta_v_preliminary_very_good(w)= delta_v_preliminary_good(X,Y,Z);
                          delta_v_preliminary_very_good(delta_v_preliminary_very_good==0)=NaN;
                          X_vect(w)=X;
                          Y_vect(w)=Y;
                          Z_vect(w)=Z;
                          
                          delta_v1_good_vect(w) = delta_v_1_good(X,Y);
                  
                          delta_v_3_good_vect(w)=delta_v_3_good(Y,Z);

                          w=w+1;
                       end
        end
    end
end

% minimizing the accurate solutions for each window formed by the inaccurate solutions
[delta_v_preliminary_min(numero),index_min]=min(delta_v_preliminary_very_good);
day_min_departure_in_window=X_vect(index_min);
day_min_arrivalSAT_in_window=Y_vect(index_min);
day_min_arrivalNEO_in_window=Z_vect(index_min);

% computing days corresponding to accurate solutions
mjd_departure_EARTH_oneday=mjd_departure_MIN_EARTH + day_min_departure_in_window*unit_good;
mjd_arrival_SAT_oneday=mjd_arrival_MIN_SAT + day_min_arrivalSAT_in_window*unit_good;
mjd_arrival_NEO_oneday=mjd_arrival_MIN_NEO + day_min_arrivalNEO_in_window*unit_good;


%% STEP 3: FLYBY
% Design of a powered flyby for each accurate solution (with min delta_v_primary) such that the spacecraft moves from the first
% to the second Lambert arc. 

mjd_departure_Earth=mjd_departure_EARTH_oneday;
mjd_arrival_SAT=mjd_arrival_SAT_oneday;
mjd_arrival_NEO=mjd_arrival_NEO_oneday;

%EARTH ORBIT in the flyby day for each region

mjd2000_departure_flyby=mjd_departure_Earth - 51544.5;
ibody_E_vect=3;
[kep_E_flyby,ksun] = uplanet(mjd2000_departure_flyby, ibody_E_vect);
[rr_0_E_flyby,vv_0_E_flyby]=par2car(kep_E_flyby(1),kep_E_flyby(2),kep_E_flyby(3),kep_E_flyby(4),kep_E_flyby(5),kep_E_flyby(6),ksun);
mu_Sun=astroConstants(4);
dth=0.1*pi/180;

%NEO ORBIT in the flyby day for each region

mjd_2000_arrival_NEO_flyby=mjd_arrival_NEO - 51544.5;
id=88;
[kep_NEO_flyby,mass_NEO,M] = ephNEO(mjd_2000_arrival_NEO_flyby,id);
[rr_0_NEO_flyby,vv_0_NEO_flyby]=par2car(kep_NEO_flyby(1),kep_NEO_flyby(2),kep_NEO_flyby(3),kep_NEO_flyby(4),kep_NEO_flyby(5),kep_NEO_flyby(6),mu_Sun);

%SATURN ORBIT in the flyby day for each region

mjd_2000_arrival_SAT_flyby=mjd_arrival_SAT - 51544.5;
ibody_Sat_vect=6;
[kep_Sat_flyby,ksun] = uplanet(mjd_2000_arrival_SAT_flyby, ibody_Sat_vect);
[rr_0_SAT_flyby,vv_0_SAT_flyby]=par2car(kep_Sat_flyby(1),kep_Sat_flyby(2),kep_Sat_flyby(3),kep_Sat_flyby(4),kep_Sat_flyby(5),kep_Sat_flyby(6),mu_Sun);

%LAMBERT TRANSFERS IN THE DAY OF FLYBY

[delta_v_1_flyby,v_final_Lambert_1_x_flyby,v_final_Lambert_1_y_flyby,v_final_Lambert_1_z_flyby] = Lambert_opt(mjd_departure_Earth,mjd_departure_Earth,mjd_arrival_SAT,mjd_arrival_SAT,unit_good,Planet_1,Planet_2,orbitType,delta_v_1_MAX);

[delta_v_3_flyby,v_initial_Lambert_3_x_flyby,v_initial_Lambert_3_y_flyby,v_initial_Lambert_3_z_flyby] = Lambert_opt_NEO(mjd_arrival_SAT,mjd_arrival_SAT,mjd_arrival_NEO,mjd_arrival_NEO,unit_good,Planet_2,NEO,orbitType,delta_v_1_MAX);

% delta_v required by the flyby
delta_v_flyby_x=(v_initial_Lambert_3_x_flyby - v_final_Lambert_1_x_flyby);
delta_v_flyby_y=(v_initial_Lambert_3_y_flyby - v_final_Lambert_1_y_flyby);
delta_v_flyby_z=(v_initial_Lambert_3_z_flyby - v_final_Lambert_1_z_flyby);

delta_v_flyby(:,numero)=[delta_v_flyby_x;delta_v_flyby_y;delta_v_flyby_z];

%Design the powered flyby around Saturn given the heliocentric velocities before and after the flyby, and Saturn's position.

V_minus_SC=[v_final_Lambert_1_x_flyby ; v_final_Lambert_1_y_flyby ; v_final_Lambert_1_z_flyby]; % vel of SC wrt SUN before flyby
V_plus_SC=[v_initial_Lambert_3_x_flyby ; v_initial_Lambert_3_y_flyby ; v_initial_Lambert_3_z_flyby]; % vel of SC wrt SUN after flyby

rr_SAT=norm(rr_0_SAT_flyby);

% compute v_minus_infinite (entry) and v_plus_infinite (exit) (wrt Saturn)

v_minus_infinite= V_minus_SC - vv_0_SAT_flyby ;
v_plus_infinite= V_plus_SC - vv_0_SAT_flyby;
delta_v=v_plus_infinite - v_minus_infinite;

v_minus_infinite_norm=norm(v_minus_infinite);
v_plus_infinite_norm=norm(v_plus_infinite);
% turn angle required by the flyby
turn_angle(numero)=acos(dot(v_minus_infinite,v_plus_infinite)/(v_plus_infinite_norm*v_minus_infinite_norm));

% in GA flyby, hyperbola is formed by 2 different hyperbola's arcs
% first arc: ENTRY ARC
mu_SAT=astroConstants(16);
e_minus_fun= @(r_P) 1 + (r_P * v_minus_infinite_norm^2)/mu_SAT;
delta_minus_fun= @(r_P) 2*asin(1./e_minus_fun(r_P));

% second arc: EXIT ARC
e_plus_fun= @(r_P) 1 + (r_P * v_plus_infinite_norm^2)/mu_SAT;
delta_plus_fun= @(r_P) 2*asin(1./e_plus_fun(r_P));

delta= @(r_P) delta_minus_fun(r_P)./2 + delta_plus_fun(r_P)./2 ; %total turn angle
fun=@(r_P) delta(r_P) - turn_angle(numero); %find the zero to fin r_P

% search for r_P: distance from the centre of Saturn to the pericentre of the flyby hyperbolas, with the constrains:
% 1: r_P > R_SAT such that the spacecraft doesn't impact the Planet during the gravity assist
% 2: r_P < R_rings_in  && r_P > R_rings_fin  such that the spacecraft doesn't hit the Rings of the Planet during the flyby
% 3: r_P < R_SOI_SAT  such that the spacecraft is inside the sphere of influence of the Planet during the flyby
% the r_P_guess has be founded using a tuning method such that the solution of fzero converges in every case 


R_SAT=astroConstants(26);
R_rings_in=R_SAT + 6600;
R_rings_fin=R_SAT + 120000;
mass_SAT=5.683 * 1e26;
mass_SUN=1.989 * 1e30;
SOI_SAT=kep_Sat_flyby(1)*(mass_SAT/mass_SUN)^(2/5);
h_atm=0;
r_P_min=R_SAT + h_atm; %minimun possible perigee for hyperbola
r_P_guess=450000;

r_P_find(numero)=fzero(fun,r_P_guess);

control(numero)=fun(r_P_find(numero));

if r_P_find(numero) < R_SAT   &&  r_P_find(numero) > SOI_SAT
    r_P_find(numero)=NaN;
else
    r_P_find(numero)=r_P_find(numero);
end

if r_P_find(numero) < R_rings_fin   &&   r_P_find(numero) > R_rings_in
    r_P_find(numero)=NaN;
else
    r_P_find(numero)=r_P_find(numero);
end

v_P_minus=sqrt(v_minus_infinite_norm^2 + (2*mu_SAT/r_P_find(numero))); %vel at pericenter at entry arc
v_P_plus=sqrt(v_plus_infinite_norm^2 + (2*mu_SAT/r_P_find(numero))); %vel at percinter at exit arc

% delta_v to give at pericentre of flyby hyperbolas
delta_v_P(numero)=abs(v_P_plus - v_P_minus);


v_minus_infinite_matr(:,numero)= V_minus_SC - vv_0_SAT_flyby ;
v_plus_infinite_matr(:,numero)= V_plus_SC - vv_0_SAT_flyby;
v_minus_infinite_norm_matr(numero)=norm(v_minus_infinite_matr(:,numero));
v_plus_infinite_norm_matr(numero)=norm(v_plus_infinite_matr(:,numero));
vv_SAT_matr(:,numero)=vv_0_SAT_flyby;
V_minus_SC_matr(:,numero)=V_minus_SC;
V_plus_SC_matr(:,numero)=V_plus_SC;
end

%% STEP 4: CONCLUSION
% Computation and minimization of the total delta_v needed for the transfer, computation of the days for departure from Earth, 
% arrival on Saturn (flyby day) and arrival on NEO. Collection of the significant parameters of the transfer.

% total and accurate delta_v needed for the whole transfer, for each region 
delta_v_TOT_good=delta_v_preliminary_min + delta_v_P;

% real min: minimum of total and accurate delta_v
[delta_v_TOT_MIN,index_VERY_MIN]=min(delta_v_TOT_good);

%find the days of delta_v_min transfer
mjd_departure_MINIMO_EARTH=mjd_departure_in + (K_vect(index_VERY_MIN)-1).*n_ES.*T_syn_ES + (I_vect(index_VERY_MIN) -1).*unit + X_vect(index_VERY_MIN)*unit_good;
mjd_arrival_MINIMO_SAT=(mjd_departure_in + n_ES.*T_syn_ES) + (L_vect(index_VERY_MIN)-1).*n_ES.*T_syn_ES + (J_vect(index_VERY_MIN) -1).*unit + Y_vect(index_VERY_MIN)*unit_good;
mjd_arrival_MINIMO_NEO=(mjd_departure_in + 2.*n_ES.*T_syn_ES) + (Q_vect(index_VERY_MIN)-1).*n_ES.*T_syn_ES + (P_vect(index_VERY_MIN) -1).*unit + Z_vect(index_VERY_MIN)*unit_good;


MJD2000_DEPARTURE=mjd_departure_MINIMO_EARTH - 51544.5;
MJD2000_FLYBY=mjd_arrival_MINIMO_SAT - 51544.5;
MJD2000_ARRIVAL=mjd_arrival_MINIMO_NEO - 51544.5;

% %% GIVEN THE BEST DAYS (DEPARTURE FROM EARTH, FLYBY AND ARRIVAL ON NEO) TO PERFORM THE MISSION, 
% % DIVIDE EACH DAY IN 24 PARTS IN ORDER TO FIND THE BEST SOLUTION UP TO 1 HOUR PRECISION.
% unit_final=unit_good/24;
% departure_final=mjd_departure_MINIMO_EARTH:unit_final:(mjd_departure_MINIMO_EARTH + 1);
% arrival_SAT_final=mjd_arrival_MINIMO_SAT:unit_final:(mjd_arrival_MINIMO_SAT + 1);
% arrival_NEO_final=mjd_arrival_MINIMO_NEO:unit_final:(mjd_arrival_MINIMO_NEO + 1);
% 
% 
% [delta_v_1_final,mjd2000_dep_opt,mjd2000_flyby_opt] = Lambert_opt_final(departure_final(1),departure_final(end),arrival_SAT_final(1),arrival_SAT_final(end),unit_final,Planet_1,Planet_2,orbitType,delta_v_1_MAX);
% 
% [delta_v_3_final,mjd2000_flyby2_opt,mjd2000_arr_opt] = Lambert_opt_NEO_final(arrival_SAT_final(1),arrival_SAT_final(end),arrival_NEO_final(1),arrival_NEO_final(end),unit_final,Planet_2,NEO,orbitType,delta_v_3_MAX);
% 







% DATAS FOR THE DAYS CHOOSEN FOR THE TRANSFER
delta_v_preliminary=delta_v_preliminary_min(index_VERY_MIN);
delta_v_P_flyby=delta_v_P(index_VERY_MIN);
turn_angle_flyby=turn_angle(index_VERY_MIN);
r_P_flyby=r_P_find(index_VERY_MIN);
delta_v_flyby_vect=delta_v_flyby(:,index_VERY_MIN);
delta_v_flyby_magnitude=norm(delta_v_flyby_vect); %note the difference between delta_v_flyby_magnitude and delta_v_P_flyby
v_minus_infinite_norm_flyby=norm(v_minus_infinite_matr(:,index_VERY_MIN));
v_plus_infinite_norm_flyby=norm(v_plus_infinite_matr(:,index_VERY_MIN));
v_minus_infinite_flyby=v_minus_infinite_matr(:,index_VERY_MIN);
v_plus_infinite_flyby=v_plus_infinite_matr(:,index_VERY_MIN);

e_minus_fun_flyby= @(r_P) 1 + (r_P * v_minus_infinite_norm_flyby^2)/mu_SAT;
delta_minus_fun_flyby= @(r_P) 2*asin(1./e_minus_fun_flyby(r_P));
e_plus_fun_flyby= @(r_P) 1 + (r_P * v_plus_infinite_norm_flyby^2)/mu_SAT;
delta_plus_fun_flyby= @(r_P) 2*asin(1./e_plus_fun_flyby(r_P));
%direction perpendicolar to flyby plane
u_flyby=cross(v_minus_infinite_flyby,v_plus_infinite_flyby)/norm(cross(v_minus_infinite_flyby,v_plus_infinite_flyby));
i_flyby=acos(u_flyby(3));
N_flyby=cross([0 ;0; 1],u_flyby);
N_flyby_norm=N_flyby/norm(N_flyby);
OM_flyby=acos(N_flyby_norm(1));
if N_flyby_norm(2) < 0
    OM_flyby=2*pi - acos(N_flyby_norm(1));
end
vv_SAT_flyby=vv_SAT_matr(:,index_VERY_MIN);
V_minus_SC_flyby=V_minus_SC_matr(:,index_VERY_MIN);
V_plus_SC_flyby=V_plus_SC_matr(:,index_VERY_MIN);
V_minus_SC_flyby_norm=norm(V_minus_SC_flyby);
V_plus_SC_flyby_norm=norm(V_plus_SC_flyby);
%leading side flyby: heliocentric velocity is reduced

%% PLOT THE ORBITS

% EARTH ORBIT AT REAL DEPARTURE DAY
mjd2000_departure_FLYBY=mjd_departure_MINIMO_EARTH - 51544.5;
[kep_E_FLYBY,ksun] = uplanet(mjd2000_departure_FLYBY, ibody_E_vect);
[rr_0_E_FLYBY,vv_0_E_FLYBY]=par2car(kep_E_FLYBY(1),kep_E_FLYBY(2),kep_E_FLYBY(3),kep_E_FLYBY(4),kep_E_FLYBY(5),kep_E_FLYBY(6),ksun);

%SATURN ORBIT AT REAL FLYBY DAY
mjd2000_arrival_SAT_FLYBY=mjd_arrival_MINIMO_SAT - 51544.5;
[kep_Sat_FLYBY,ksun] = uplanet(mjd2000_arrival_SAT_FLYBY, ibody_Sat_vect);
[rr_0_SAT_FLYBY,vv_0_SAT_FLYBY]=par2car(kep_Sat_FLYBY(1),kep_Sat_FLYBY(2),kep_Sat_FLYBY(3),kep_Sat_FLYBY(4),kep_Sat_FLYBY(5),kep_Sat_FLYBY(6),mu_Sun);

%NEO ORBIT AT REAL ARRIVAL DAY
mjd2000_arrival_NEO_FLYBY=mjd_arrival_MINIMO_NEO - 51544.5;
[kep_NEO_FLYBY,mass_NEO,M] = ephNEO(mjd2000_arrival_NEO_FLYBY,id);
[rr_0_NEO_FLYBY,vv_0_NEO_FLYBY]=par2car(kep_NEO_FLYBY(1),kep_NEO_FLYBY(2),kep_NEO_FLYBY(3),kep_NEO_FLYBY(4),kep_NEO_FLYBY(5),kep_NEO_FLYBY(6),mu_Sun);

AU=astroConstants(2);
% plot Earth orbit
figure
plotOrbit(kep_E_FLYBY(1),kep_E_FLYBY(2),kep_E_FLYBY(3),kep_E_FLYBY(4),kep_E_FLYBY(5),kep_E_FLYBY(6),kep_E_FLYBY(6)+2*pi,dth,mu_Sun);
hold on
%plot Earth
plot3(rr_0_E_FLYBY(1),rr_0_E_FLYBY(2),rr_0_E_FLYBY(3),'bo','LineWidth',3)
hold on
%plot NEO orbit
plotOrbit(kep_NEO_FLYBY(1),kep_NEO_FLYBY(2),kep_NEO_FLYBY(3),kep_NEO_FLYBY(4),kep_NEO_FLYBY(5),kep_NEO_FLYBY(6),kep_NEO_FLYBY(6)+2*pi,dth,mu_Sun);
hold on
%plot NEO
hold on
plot3(rr_0_NEO_FLYBY(1),rr_0_NEO_FLYBY(2),rr_0_NEO_FLYBY(3),'ko','LineWidth',3)
%plot Saturn orbit
plotOrbit(kep_Sat_FLYBY(1),kep_Sat_FLYBY(2),kep_Sat_FLYBY(3),kep_Sat_FLYBY(4),kep_Sat_FLYBY(5),kep_Sat_FLYBY(6),kep_Sat_FLYBY(6)+2*pi,dth,mu_Sun);
%plot Saturn
hold on
plot3(rr_0_SAT_FLYBY(1),rr_0_SAT_FLYBY(2),rr_0_SAT_FLYBY(3),'ro','LineWidth',3)
%plot the Sun
hold on
plot3(0,0,0,'yo','LineWidth',5)

%% PLOT THE TRANSFERS

Nrev=0; 
optionsLMR=1;
orbitType=0;

% Transfer between Earth and Saturn
TOF_T1=(mjd_arrival_MINIMO_SAT-mjd_departure_MINIMO_EARTH)*24*60*60;
[a_T1,p_T1,e_T1,ERROR,v_i_T1,v_f_T1,TPAR_T1,THETA_T1] = lambertMR(rr_0_E_FLYBY,rr_0_SAT_FLYBY,TOF_T1,mu_Sun,orbitType,Nrev,optionsLMR);

T_T1=(2*pi)*sqrt(a_T1^3/ksun);
tspan_T1=linspace(0,TOF_T1,1000);

%options
options = odeset('RelTol',1e-13,'AbsTol',1e-14);

[time_T1, state_T1 ] = ode45(@(t,s) twobody_problem_ode(t,s,ksun), tspan_T1, [rr_0_E_FLYBY v_i_T1'],options);

%plot the transfer orbit

hold on
plot3(state_T1(:,1), state_T1(:,2), state_T1(:,3) , 'r', 'LineWidth', 2);


% Transfer between Saturn and NEO
TOF_T2=(mjd_arrival_MINIMO_NEO-mjd_arrival_MINIMO_SAT)*24*60*60;
[a_T2,p_T2,e_T2,ERROR,v_i_T2,v_f_T2,TPAR_T2,THETA_T2] = lambertMR(rr_0_SAT_FLYBY,rr_0_NEO_FLYBY,TOF_T2,mu_Sun,orbitType,Nrev,optionsLMR);

T_T2=(2*pi)*sqrt(a_T2^3/ksun);
tspan_T2=linspace(0,TOF_T2,1000);

%options
options = odeset('RelTol',1e-13,'AbsTol',1e-14);

[time_T2, state_T2 ] = ode45(@(t,s) twobody_problem_ode(t,s,ksun), tspan_T2, [rr_0_SAT_FLYBY v_i_T2'],options);

%plot the transfer orbit

hold on
plot3(state_T2(:,1), state_T2(:,2), state_T2(:,3) , 'k', 'LineWidth', 2);
legend('Earth orbit','Earth','NEO orbit','NEO','Saturn orbit','Saturn','Sun','1st Leg','2nd Leg')
title 'Earth-Saturn-NEO transfer'

%% PLOT THE HYPERBOLA ARCS OF THE FLYBY

e_minus=e_minus_fun_flyby(r_P_flyby);
e_plus=e_plus_fun_flyby(r_P_flyby);
a_minus=-mu_SAT/v_minus_infinite_norm_flyby^2;
a_plus=-mu_SAT/v_plus_infinite_norm_flyby^2;
delta_minus=delta_minus_fun_flyby(r_P_flyby);
delta_plus=delta_plus_fun_flyby(r_P_flyby);


% PLOT HYPERBOLAS IN PLANETOCENTRIC FRAME
% INCOMING HYPERBOLA
beta_minus=(pi - delta_minus)/2;
i_minus=i_flyby;
om_minus=beta_minus;
OM_minus=OM_flyby;
theta_0_in=-0.5*pi;
theta_f_in=0;
dth=0.5*pi/180;


%PLOT SATURN
opts_example2.RefPlane = 'ecliptic';
figure;
planet3D('Saturn',opts_example2);
grid on;

%PLOT INCOMING HYPERBOLA
if(theta_f_in<theta_0_in)
   theta_f_in=theta_0_in+2*pi; 
end
[rr_in] = par2car(a_minus,e_minus,i_minus,OM_minus,om_minus, theta_0_in,mu_SAT);
Posizioni_in=1000*[rr_in]; %from [Km] to [m]
for th=theta_0_in:dth:theta_f_in
[rr_in] = par2car(a_minus,e_minus,i_minus,OM_minus,om_minus, th,mu_SAT);
Posizioni_in=[Posizioni_in,1000*rr_in];
end
hold on
plot3(Posizioni_in(1,:),Posizioni_in(2,:),Posizioni_in(3,:),'r','LineWidth',2)
axis equal;

% OUTCOMING HYPERBOLA
beta_plus=(pi - delta_plus)/2;
i_plus=i_flyby;
om_plus=beta_plus;
OM_plus=OM_flyby;
theta_0_out=0;
theta_f_out=0.5*pi;
dth=0.5*pi/180;

%PLOT OUTCOMING HYPERBOLA
if(theta_f_out<theta_0_out)
   theta_f_out=theta_0_out+2*pi; 
end
[rr_out] = par2car(a_plus,e_plus,i_plus,OM_plus,om_plus, theta_0_out,mu_SAT);
Posizioni_out=1000*[rr_out]; %from [Km] to [m]
for th=theta_0_out:dth:theta_f_out
[rr_out] = par2car(a_plus,e_plus,i_plus,OM_plus,om_plus, th,mu_SAT);
Posizioni_out=[Posizioni_out,1000*rr_out];
end
hold on
plot3(Posizioni_out(1,:),Posizioni_out(2,:),Posizioni_out(3,:),'k','LineWidth',2)
hold on
quiver3(0,0,0,vv_0_SAT_FLYBY(1)*2e7,vv_0_SAT_FLYBY(2)*2e7,vv_0_SAT_FLYBY(3)*2e7,'LineWidth',1.5)

xlabel 'x [m]'
ylabel 'y [m]'
zlabel 'z [m]'

legend('','incoming hyperbola','outcoming hyperbola','Saturn velocity')
title 'flyby around Saturn'

%% STEP 5: TIME DURATION OF THE FLYBY
% In order to solve the problem in the first place,the patched conics method has been used, wich implies that the sphere of influence of the
% planet is considered to be infinitesimal with reference to the heliocentric legs; therefore the time duration of the flyby is considered to be zero (impulsive manouvre).
% Instead, when the problem is treated from the planetary point of view (as has be done in order to design the powered gravity assist),
% the sphere of influence of the planet can be considerate to be finite (and can be computed).
% In this framework, the time duration of the flyby can becomputed as the time passing from the entry to the exit in/from the sphere of influence of the Planet.


theta_inf_minus=acos(-1/e_minus);
theta_inf_plus=acos(-1/e_plus);
% INCOMING HYPERBOLA: time passing from the entry inside SOI_SAT, until the spacecraft reaches the pericentre of the hyperbola
p_minus=a_minus*(1 - e_minus^2);
theta_minus_SOI=acos((p_minus - SOI_SAT)/(SOI_SAT*e_minus));
F_minus=2*atanh(sqrt((e_minus-1)/(1+e_minus)) * tan((theta_minus_SOI)/2));

delta_t_minus=sqrt(-a_minus^3/mu_SAT) * (e_minus * sinh(F_minus) - F_minus);


% OUTCOMING HYPERBOLA: time passing from the passage at the pericentre of the Hyperbola, until the exit from the SOI
p_plus=a_plus*(1 - e_plus^2);
theta_plus_SOI=acos((p_plus - SOI_SAT)/(SOI_SAT*e_plus));
F_plus=2*atanh(sqrt((e_plus-1)/(1+e_plus)) * tan(( theta_plus_SOI)/2));

delta_t_plus=sqrt(-a_plus^3/mu_SAT) * (e_plus * sinh(F_plus) - F_plus);

TOF_flyby=delta_t_minus + delta_t_plus;
TOF_flyby_days=TOF_flyby/(60*60*24);








 case 'DIRECT TRANSFER'
%% DIRECT TRANSFER BETWEEN EARTH AND NEO


clc
clear
close all


%EARTH ORBIT
mjd_departure_in=63139;
mjd2000_departure=mjd_departure_in - 51544.5;
ibody_E_vect=3;
[kep_E_in,ksun] = uplanet(mjd2000_departure, ibody_E_vect);
mu_Sun=astroConstants(4);
%NEO ORBIT
id=88;
[kep_NEO_in,mass_NEO,M] = ephNEO(mjd2000_departure,id);
%SATURN ORBIT
ibody_Sat_vect=6;
[kep_Sat_in,ksun] = uplanet(mjd2000_departure, ibody_Sat_vect);

Planet_1=3;

NEO=88;
orbitType=0;
delta_v_1_MAX=75;
delta_v_3_MAX=75;
mjd_departure_in=63139;
% 1 synodic period between Earth and NEO for departure window initial guess:
n_ENEO=1;
T_E=2*pi*sqrt(kep_E_in(1).^3/mu_Sun);
T_NEO=2*pi*sqrt(kep_NEO_in(1).^3/mu_Sun);
T_syn_ENEO=((T_NEO*T_E)/abs(T_NEO-T_E))/(60*60*24);
mjd_departure_fin= mjd_departure_in + n_ENEO*T_syn_ENEO;


% from first departure to last arrival: 75737 is 28/03/66
MAX_time=75737 - 63139;

N_1=floor(MAX_time/(n_ENEO*T_syn_ENEO));
unit=10;

%l=changing of arrival window on Saturn
 %k=changing of departing window from Earth
for k=1:N_1-1  %DEPARTURE WINDOWS FROM EARTH

    mjd_departure_in_vect(k)=mjd_departure_in + (k-1)*n_ENEO*T_syn_ENEO;
    mjd_departure_fin_vect(k)= mjd_departure_in_vect(k) + n_ENEO*T_syn_ENEO;

       for l=k:N_1-1  %ARRIVAL WINDOWS ON NEO

    mjd_arrival_in_vect(l)=mjd_departure_fin + (l-1)*n_ENEO*T_syn_ENEO;
    mjd_arrival_fin_vect(l)= mjd_arrival_in_vect(l) + n_ENEO*T_syn_ENEO;

       %ARRIVING ON NEO IN EVERY ARRIVAL WINDOW DEPARTING FROM EARTH AT EVERY DEPARTURE WINDOW
    [delta_v_tot(:,:,k,l),delta_v_min(k,l),mjd_departure_optimal(:,:,k,l),mjd_arrival_optimal(:,:,k,l)] = Lambert_optimal_NEO(mjd_departure_in_vect(k),mjd_departure_fin_vect(k),mjd_arrival_in_vect(l),mjd_arrival_fin_vect(l),unit,Planet_1,NEO,orbitType,delta_v_1_MAX);

   end
end



z=1;
for I=1:size(delta_v_tot,1) % departing day from Earth
    for J=1:size(delta_v_tot,2) % arrival day on NEO 
        for K=1:size(delta_v_tot,3) % departing window from Earth
            for L=K:size(delta_v_tot,4) % arrival window on NEO 

                 if  delta_v_tot(I,J,K,L)< 12
                           
                          delta_v_TOT_good(z)= delta_v_tot(I,J,K,L);
                          I_vect(z)=I;
                          J_vect(z)=J;
                          K_vect(z)=K;
                          L_vect(z)=L;

                          z=z+1;
                 end
  
            end
        end
    end
end

[delta_v_TOT_preliminar_MIN,index_MIN]=min(delta_v_TOT_good);
% FOR CYCLE FOR ANY GOOD SOLUTION WICH RESPECT THE CONSTRAIN
% delta_v_TOT_good_preliminar< 16.5


mjd_departure_MIN_EARTH=mjd_departure_in + (K_vect(index_MIN)-1).*n_ENEO.*T_syn_ENEO + (I_vect(index_MIN) -1).*unit;
mjd_arrival_MIN_NEO=(mjd_departure_in + n_ENEO.*T_syn_ENEO) + (L_vect(index_MIN)-1).*n_ENEO.*T_syn_ENEO + (J_vect(index_MIN) -1).*unit;

 unit_good=1;
 
% good solution discretized day by day
% windows 
mjd_departure_EARTH_WINDOW=mjd_departure_MIN_EARTH:unit_good:(mjd_departure_MIN_EARTH+unit);
mjd_arrival_NEO_WINDOW=mjd_arrival_MIN_NEO:unit_good:(mjd_arrival_MIN_NEO+unit);


[delta_v_tot_final,delta_v_min_final,mjd_departure_optimal_final,mjd_arrival_optimal_final,delta_v_1,delta_v_2] = Lambert_optimal_NEO(mjd_departure_EARTH_WINDOW(1),mjd_departure_EARTH_WINDOW(end),mjd_arrival_NEO_WINDOW(1),mjd_arrival_NEO_WINDOW(end),unit_good,Planet_1,NEO,orbitType,delta_v_1_MAX);


%% PLOT ORBITS AND TRANSFER
%% PLOT THE ORBITS
dth=0.5*pi/180;
% EARTH ORBIT AT DEPARTURE DAY
mjd2000_departure_opt=mjd_departure_optimal_final - 51544.5;
[kep_E,ksun] = uplanet(mjd2000_departure_opt, ibody_E_vect);
[rr_0_E,vv_0_E]=par2car(kep_E(1),kep_E(2),kep_E(3),kep_E(4),kep_E(5),kep_E(6),ksun);

%NEO ORBIT AT ARRIVAL DAY
mjd2000_arrival_NEO_opt=mjd_arrival_optimal_final - 51544.5;
[kep_NEO,mass_NEO,M] = ephNEO(mjd2000_arrival_NEO_opt,id);
[rr_0_NEO,vv_0_NEO]=par2car(kep_NEO(1),kep_NEO(2),kep_NEO(3),kep_NEO(4),kep_NEO(5),kep_NEO(6),mu_Sun);


%plot Earth's orbit
plotOrbit(kep_E(1),kep_E(2),kep_E(3),kep_E(4),kep_E(5),kep_E(6),kep_E(6)+2*pi,dth,mu_Sun);
%plot NEO's orbit
hold on
plotOrbit(kep_NEO(1),kep_NEO(2),kep_NEO(3),kep_NEO(4),kep_NEO(5),kep_NEO(6),kep_NEO(6)+2*pi,dth,mu_Sun);
%plot the Sun
hold on
plot3(0,0,0,'yo','LineWidth',5)
%plot Earth
plot3(rr_0_E(1),rr_0_E(2),rr_0_E(3),'bo','LineWidth',3)
%plot NEO
hold on
plot3(rr_0_NEO(1),rr_0_NEO(2),rr_0_NEO(3),'ko','LineWidth',3)

title 'direct transfer Earth-NEO'

%% PLOT THE TRANSFER

%	orbitType[1]    Logical variable defining whether transfer is
%                       0: direct transfer from R1 to R2 (counterclockwise)
Nrev=0; %number of revolutions
optionsLMR=1;%warnings are displayed only when the algorithm does not converge
orbitType=0;

% Transfer between Earth and NEO
TOF_T1=([mjd_arrival_optimal_final-mjd_departure_optimal_final])*24*60*60;
[a_T1,p_T1,e_T1,ERROR,v_i_T1,v_f_T1,TPAR_T1,THETA_T1] = lambertMR(rr_0_E,rr_0_NEO,TOF_T1,mu_Sun,orbitType,Nrev,optionsLMR);

T_T1=(2*pi)*sqrt(a_T1^3/ksun);
tspan_T1=linspace(0,TOF_T1,1000);

%options
options = odeset('RelTol',1e-13,'AbsTol',1e-14);

[time_T1, state_T1 ] = ode45(@(t,s) twobody_problem_ode(t,s,ksun), tspan_T1, [rr_0_E,v_i_T1'],options);

%plot the transfer orbit

hold on
plot3(state_T1(:,1), state_T1(:,2), state_T1(:,3) , 'k', 'LineWidth', 1.5);

legend('Earth orbit','NEO orbit','Sun','Earth','NEO','Lambert leg')

end