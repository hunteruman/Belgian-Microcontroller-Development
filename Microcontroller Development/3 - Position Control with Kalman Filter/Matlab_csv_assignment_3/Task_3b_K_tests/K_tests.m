clear all
close all
clear global
clc

%% Constants
fs = 100; % sampling frequency (control function - 100Hz)
Ts = 1/fs;

%% Data input and plot
pref_col = 13;
dist_col = 9;
datafile1 = 'PositionStepK05.csv'; % K=0.5
datafile2 = 'PositionStepK1.csv'; % K=1
datafile3 = 'PositionStepK2.csv'; % K=2
datafile4 = 'PositionStepK4.csv'; % K=4
Metafile1 = csvread(datafile1); % M_t and E_t
Metafile2 = csvread(datafile2); % M_t and E_t
Metafile3 = csvread(datafile3); % M_t and E_t
Metafile4 = csvread(datafile4); % M_t and E_t
% -------------------------------------------
N = length(Metafile4); % = 2179 (K=4), adapt according to the smallest file size
pref1_t = Metafile1(1:N,pref_col); % Step reference
dist1_t = -Metafile1(1:N,dist_col); % Measured response
pref2_t = Metafile2(1:N,pref_col);
dist2_t = -Metafile2(1:N,dist_col);
pref3_t = Metafile3(1:N,pref_col);
dist3_t = -Metafile3(1:N,dist_col);
pref4_t = Metafile4(1:N,pref_col);
dist4_t = -Metafile4(1:N,dist_col);
t = [0:N-1]'*Ts; % time elapsed

% Overall plot
figure
subplot(2,2,1),plot(t, [pref1_t, dist1_t], 'LineWidth', 1)
title('Position step reference vs measured response for K=0.5')
grid on
axis tight
xlabel('t [s]')
ylabel('Pos [m]')
legend('Measured', 'Estimate')

subplot(2,2,2),plot(t, [pref2_t, dist2_t], 'LineWidth', 1)
title('Position step reference vs measured response for K=1')
grid on
axis tight
xlabel('t [s]')
ylabel('Pos [m]')
legend('Measured', 'Estimate')

subplot(2,2,3),plot(t, [pref3_t, dist3_t], 'LineWidth', 1)
title('Position step reference vs measured response for K=2')
grid on
axis tight
xlabel('t [s]')
ylabel('Pos [m]')
legend('Measured', 'Estimate')

subplot(2,2,4),plot(t, [pref4_t, dist4_t], 'LineWidth', 1)
title('Position step reference vs measured response for K=4')
grid on
axis tight
xlabel('t [s]')
ylabel('Pos [m]')
legend('Measured', 'Estimate')

% Detailed plot

N2 = 100;
pref1_t = Metafile1(1:N2,pref_col); % Step reference
dist1_t = -Metafile1(1:N2,dist_col); % Measured response
pref2_t = Metafile2(1:N2,pref_col);
dist2_t = -Metafile2(1:N2,dist_col);
pref3_t = Metafile3(1:N2,pref_col);
dist3_t = -Metafile3(1:N2,dist_col);
pref4_t = Metafile4(1:N2,pref_col);
dist4_t = -Metafile4(1:N2,dist_col);
t = [0:N2-1]'*Ts; % time elapsed

% Overall plot
figure
subplot(2,2,1),plot(t, [pref1_t, dist1_t], 'LineWidth', 1)
title('Measured vs estimate position plot Q/R=0.001')
grid on
axis tight
xlabel('t [s]')
ylabel('Pos [m]')
legend('Measured', 'Estimate')

subplot(2,2,2),plot(t, [pref2_t, dist2_t], 'LineWidth', 1)
title('Measured vs estimate position plot Q/R=0.01')
grid on
axis tight
xlabel('t [s]')
ylabel('Pos [m]')
legend('Measured', 'Estimate')

subplot(2,2,3),plot(t, [pref3_t, dist3_t], 'LineWidth', 1)
title('Measured vs estimate position plot Q/R=0.1')
grid on
axis tight
xlabel('t [s]')
ylabel('Pos [m]')
legend('Measured', 'Estimate')

subplot(2,2,4),plot(t, [pref4_t, dist4_t], 'LineWidth', 1)
title('Measured vs estimate position plot Q/R=1')
grid on
axis tight
xlabel('t [s]')
ylabel('Pos [m]')
legend('Measured', 'Estimate')
%display ('Press any key to continue '); pause;