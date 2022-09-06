clear all; 
close all; 
format short e; 

% Simulation parameters :
% -------------------------

f = 3.6e9; % frequency [GHz]
c = 3e8; % light velocity [m]
d = c/f; %distance between two antennas [m]

step = 1000-1; % set the number of points for the graphical representation
theta = [0:pi/step:pi];
gamma = pi*sin(theta);

M=input('number of antenna : M >> '); 

% M*1 steering vector
S = zeros(M,length(theta));
for m = 0:1:(M-1)
    S(m+1,:)=exp(-j*m*gamma);
end;

% Generation of the Vandermolde matrix
V = zeros(M,M);
for a = 0:1:(M-1)
    for b = 0:1:(M-1)
        V(a+1,b+1)=exp(j*2*pi*a*b/M);
    end;
end;

% Computation of the gain
gain_m = cell(M,1);
y_m = cell(M,1);
gain_m_dB = cell(M,1);
for k =1:1:M
    gain_m{k} = (abs((1/sqrt(M)*V(k,:)*S))).^2;
    y_m{k} = circshift(gain_m{k},length(gain_m{k})/2);
    gain_m_dB{k} = 10*log10(y_m{k});
end;



% g vs linear theta
figure(1);
x = [-90:(180/step):90];
hold on;
for k=1:1:M
    plot(x,y_m{k},'-','LineWidth',1);
end;
hold off;
title({"ULA with M = "+ num2str(M) + "  Vandermolde matrix, " + num2str(M) + " orthogonal beams"});
axis([-90 90 0 M]);
xlabel('angle [deg]');
ylabel('gain');
grid;

% g[dB] vs linear theta
figure(2);
hold on;
for k=1:1:M
    x = [-90:(180/step):90];
    plot(x,gain_m_dB{k},'-','LineWidth',1);
end;
hold off;
title({"ULA with M = "+ num2str(M) + "  Vandermolde matrix, " + num2str(M) + " orthogonal beams"});
axis([-90 90 -40 2*M]);
xlabel('angle [deg]');
ylabel('gain [dB]');
grid;

% g vs polar theta
figure(3);
theta = linspace(0,pi,step+1);
polarplot(theta, gain_m{1},'b-',-theta,gain_m{1},'b-', 'LineWidth', 2);
for k=2:1:M
hold on;
polarplot(theta, gain_m{k},-theta, gain_m{k},'-', 'LineWidth', 2);
hold off;
end;
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'counterclockwise';
ax.ThetaLim = [-90 90];
title({"ULA with M = "+ num2str(M) + "  Vandermolde matrix, " + num2str(M) + " orthogonal beams"});

