clear all; 
close all; 
format short e; 

% Simulation parameters :
% -------------------------

f = 3.6e9; % frequency [GHz]
c = 3e8; % light velocity [m]
d = c/f; %distance between two antennas [m]

M = [4 8 16 32 64];

step = 1000-1; % set the number of points for the graphical representation
theta = [0:pi/step:pi];
gamma = pi*sin(theta);

% M*1 steering vector
Sm = cell(length(M),1);
for k = 1:1:length(M)
    S = zeros(M(k),length(theta));
    for m = 0:1:(M(k)-1)
        S(m+1,:)=exp(-j*m*gamma);
    end;
    Sm{k} = S; 
end;

% Generation of 1*M beamformer vector
Bm = cell(length(M),1);
for k=1:1:length(M)
    one = zeros(1,M(k));
    one(1:end)=1;
    B = (1/sqrt(M(k)))*one;
    Bm{k} = B;
end;

% Computation of the gain
gain_m = cell(length(M),1);
for k =1:1:length(M)
    gain = (abs(Bm{k}*Sm{k})).^2;
    y = circshift(gain(1,:),length(gain(1,:))/2);
    gain_m{k} = 10*log10(y/max(y));
end;

% Computation of the main lobe width
coordinate = [];
x = [-90:(180/step):90];
x_coordinate = [];
width = [];

for p=1:1:length(M)
    local_min = islocalmin(gain_m{p});
    for k=1:1:length(local_min)
        if local_min(k)==1
            coordinate(end+1) = k;
            x_coordinate(end+1) = x(k);
        end
    end;
    w = min(abs(x_coordinate))*2; % because it is symmetrical
    coordinate=[];
    x_coordinate=[];
    width(end+1)=w;
end;



% normalized g[dB] vs linear theta
figure(1);
hold on;
x = [-90:(180/step):90];
plot(x,gain_m{1},'b','LineWidth',1.5);
plot(x,gain_m{2},'r','LineWidth',1.5);
plot(x,gain_m{3},'y','LineWidth',1.5);
plot(x,gain_m{4},'m','LineWidth',1.5);
plot(x,gain_m{5},'g','LineWidth',1.5);
hold off;
title("Normalized gain g/max(g) [dB]");
axis([-90 90 -30 max(gain_m{1})]);
xlabel('angle [deg]');
ylabel('gain [dB]');
legend('4','8','16','32','64');
grid;

figure(2);
plot(M,width,'b-o','LineWidth',1.5);
axis([0 64 0 max(width)]);
xlabel('main lobe width [deg]');
ylabel('number of antennas M');
grid;