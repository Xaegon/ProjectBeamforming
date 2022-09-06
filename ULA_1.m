clear all; 
close all; 
format short e; 

% Simulation parameters :
% -------------------------

f = 3.6e9; % frequency [GHz]
c = 3e8; % light velocity [m]
d = c/f; %distance between two antennas [m]

M=input('number of antenna : M >> '); 

step = 1000-1; % set the number of points for the graphical representation
theta = [0:pi/step:pi];
gamma = pi*sin(theta);

% M*1 steering vector
S = zeros(M,length(theta));
for m = 0:1:(M-1)
    S(m+1,:)=exp(-j*m*gamma);
end;

% Generation of 1*M beamformer vector
one = zeros(1,M);
one(1:end)=1;
B = (1/sqrt(M))*one;

% Computation of the gain
gain = (abs(B*S)).^2;
y = circshift(gain(1,:),length(gain(1,:))/2);



% g vs linear theta
figure(1);
x = [-90:(180/step):90];
plot(x,y,'r','LineWidth',1.5);
##title({"ULA with M = "+ num2str(M) + "  radiation pattern"});
axis([-90 90 0 max(gain(1,:))]);
xlabel('angle [deg]')
ylabel('gain')
grid;

% g[dB] vs linear theta
figure(2);
x = [-90:(180/step):90];
plot(x,10*log10(y),'r','LineWidth',1.5);
##title({"ULA with M = "+ num2str(M) + "  radiation pattern"});
axis([-90 90 -40 max(10*log(y))]);
xlabel('angle [deg]')
ylabel('gain [dB]')
grid;

% g vs polar theta
figure(3);
polarplot(theta, gain(1,:), 'r',-theta, gain(1,:), 'r', 'LineWidth', 2);
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'counterclockwise';
ax.ThetaLim = [-90 90];
##title({"ULA with M = "+ num2str(M) + "  radiation pattern"});

