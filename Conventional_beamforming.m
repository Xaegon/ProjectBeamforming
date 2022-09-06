clear all; 
close all; 
format short e; 

% Simulation parameters :
% -------------------------

f = 3.6e9; % frequency [GHz]
c = 3e8; % light velocity [m]
d = c/f; %distance between two antennas [m]

step = 1000-1; % set the number of points for the graphical representation
theta = [-pi/2:pi/step:pi/2];
gamma = pi*sin(theta);

M=input('number of antenna : M >> '); 
theta0=input('given angle in degree: theta0 >> '); 
gamma0 = pi*sin(theta0*pi/180);

% 1*M beamformer vector
b = [];
for m = 0:1:(M-1)
    b(end+1)=exp(j*m*gamma0);
end;
B = (1/sqrt(M))*b ;

% M*1 steering vector
S = zeros(M,length(theta));
for m = 0:1:(M-1)
    S(m+1,:)=exp(-j*m*gamma);
end;

% Computation of the gain
gain = (abs(B*S)).^2;


% For polar chart
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = [0:pi/step:pi];
gamma = pi*sin(theta);
gamma0 = pi*sin(theta0*pi/180);

% 1*M beamformer vector
b = [];
for m = 0:1:(M-1)
    b(end+1)=exp(j*m*gamma0);
end;
B = (1/sqrt(M))*b ;

% M*1 steering vector
S = zeros(M,length(theta));
for m = 0:1:(M-1)
    S(m+1,:)=exp(-j*m*gamma);
end;

% Computation of the gain
gain_p = (abs(B*S)).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% g vs linear theta
figure(1);
x = [-90:(180/step):90];
plot(x,gain,'r','LineWidth',1.5);
title({"ULA with M = "+ num2str(M) + "  conventional beamforming \theta_0 = " + num2str(theta0) + ""});
axis([-90 90 0 max(gain(1,:))]);
xlabel('angle [deg]')
ylabel('gain')
grid;

% g[dB] vs linear theta
figure(2);
x = [-90:(180/step):90];
plot(x,10*log10(gain),'r','LineWidth',1.5);
title({"ULA with M = "+ num2str(M) + "  conventional beamforming \theta_0 = " + num2str(theta0) + ""});
axis([-90 90 -40 max(10*log(gain))]);
xlabel('angle [deg]')
ylabel('gain [dB]')
grid;

% g vs polar theta
figure(3);
polarplot([0:pi/step:pi], gain_p(1,:), 'r', 'LineWidth', 2);
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'counterclockwise';
ax.ThetaLim = [-90 90];
title({"ULA with M = "+ num2str(M) + "  conventional beamforming \theta_0 = " + num2str(theta0) + ""});