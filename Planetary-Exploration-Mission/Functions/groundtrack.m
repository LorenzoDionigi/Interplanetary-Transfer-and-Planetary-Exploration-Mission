function  [lambda,phi] =groundtrack(t,Y,omega_e)

% AUTHOR: 
 % Lorenzo Dionigi

% SUPERVISOR:
 % Prof. Camilla Colombo


delta= zeros(1,length(t)); alpha=delta;
lambda=NaN(1,length(t)); phi=lambda;
theta_G0=0;
theta_G = theta_G0 + t*omega_e;%[rad/s]
for ii = 1:length(t)
    r = norm(Y(ii,1:3));
    delta(ii)= asin(Y(ii,3)/r);
    phi(ii)=delta(ii)*180/pi;
    if Y(ii,2)/r <= 0
        alpha(ii)= 2*pi - acos ( Y(ii,1)/( r*cos(delta(ii)) ) );
    else 
        alpha(ii) = acos ( Y(ii,1)/( r*cos(delta(ii)) ) );
    end
    lambda(ii)= alpha(ii)-theta_G(ii);
    lambda(ii) = rad2deg(lambda(ii));
    lambda(ii) = wrapTo180(lambda(ii));
end

figure()
hold on
plot(lambda,phi, 'y.', 'LineStyle','none')
plot(lambda(1), phi(1), 'r o', MarkerSize=10, MarkerFaceColor='r') %plot initial point
plot(lambda(end), phi(end), 'g s ',MarkerSize=10, MarkerFaceColor='g') %plot final point

axis image;
grid on;
x_plot = [-182 182]; y_plot = [-92 92];
h=image(x_plot,-y_plot,imread('Earth_2d.jpg'));
uistack(h,'bottom')
legend('Nominal ground track perturbed - 50 period','Start','End')

% axes('position', [0.2 0.6 0.3 0.3])
% hold on
% h=image(x_plot,-y_plot,imread('Earth_2d.jpg'));
% uistack(h,'bottom')
% plot(lambda,phi, 'y.', 'LineStyle','none')
% plot(lambda(1), phi(1), 'r o', MarkerSize=10, MarkerFaceColor='r') %plot initial point
% plot(lambda(end), phi(end), 'g s ',MarkerSize=10, MarkerFaceColor='g') %plot final point
% axis([-10, 10, -10, 10])
