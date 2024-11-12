function plotOrbit(a,e,i,Om,om,th0,thf,dth,mu)

%      3D orbit plot
%
%INPUT:
 % a                major semiaxis  [Km]
 % e                eccentricity [-]
 % i                inclination [rad]
 % OM               RAAN [rad]
 % om               pericenter anomaly [rad]
 % th0              initial thrue anomaly [rad]
 % thf              final thrue anomaly [rad]
 % dth              thrue anomaly step [rad]
 % mu               planet gravitational constant   [Km^3/s^2]
 
 % CONTRIBUTORS
  % Francesco Paolo Vacca
  % Giuseppe Antonio Zito

if(thf<th0)
   thf=th0+2*pi; 
end
[rr] = par2car(a,e,i,Om,om, th0,mu);
Posizioni=[rr];
for th=th0:dth:thf
[rr] = par2car(a,e,i,Om,om, th,mu);
Posizioni=[Posizioni,rr];
end
% Terra_3D(6378)
plot3(Posizioni(1,:),Posizioni(2,:),Posizioni(3,:),'LineWidth',2)
hold on
%plot3(Posizioni(1,1),Posizioni(2,1),Posizioni(3,1),'o','MarkerEdgeColor','red')
axis equal;
grid on
xlabel 'x [Km]'
ylabel 'y [Km]'
zlabel 'z [Km]'