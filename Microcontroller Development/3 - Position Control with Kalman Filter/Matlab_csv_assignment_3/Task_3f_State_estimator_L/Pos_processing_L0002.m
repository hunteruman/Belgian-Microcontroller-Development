clear all
close all
clear global
clc

%% Constants
fs = 100; % sampling frequency (control function - 100Hz)
Ts = 1/fs;

%% Data input and plot
M_col = 9;
E_col = 10;
datafile1 = 'Badinit0001.csv'; % csvread werkt niet echt met matrices
datafile2 = 'Badinit001.csv';
datafile3 = 'PolePlacement0001.csv';
datafile4 = 'PolePlacement001.csv';
Metafile1 = csvread(datafile1); % M_t and E_t
Metafile2 = csvread(datafile2); % M_t and E_t
Metafile3 = csvread(datafile3); % M_t and E_t
Metafile4 = csvread(datafile4); % M_t and E_t
% -------------------------------------------
N = 2400; % adapted to 24s of singal
M1_t = Metafile1(1:N,M_col); % Measured
E1_t = -Metafile1(1:N,E_col); % Estimate
M2_t = Metafile2(1:N,M_col);
E2_t = -Metafile2(1:N,E_col);
M3_t = Metafile3(1:N,M_col);
E3_t = -Metafile3(1:N,E_col);
M4_t = Metafile4(1:N,M_col);
E4_t = -Metafile4(1:N,E_col);
t = [0:N-1]'*Ts; % time elapsed

% Overall plot
figure
subplot(2,2,1),plot(t, [M1_t, E1_t], 'LineWidth', 1)
title('Measured vs estimate position plot Q/R=0.001 original')
grid on
axis tight
xlabel('t [s]')
ylabel('Pos [m]')
legend('Measured', 'Estimate')

subplot(2,2,2),plot(t, [M2_t, E2_t], 'LineWidth', 1)
title('Measured vs estimate position plot Q/R=0.01 original')
grid on
axis tight
xlabel('t [s]')
ylabel('Pos [m]')
legend('Measured', 'Estimate')

subplot(2,2,3),plot(t, [M3_t, E3_t], 'LineWidth', 1)
title('Measured vs estimate position plot Q/R=0.001 L=-0.002')
grid on
axis tight
xlabel('t [s]')
ylabel('Pos [m]')
legend('Measured', 'Estimate')

subplot(2,2,4),plot(t, [M4_t, E4_t], 'LineWidth', 1)
title('Measured vs estimate position plot Q/R=0.01 L=-0.002')
grid on
axis tight
xlabel('t [s]')
ylabel('Pos [m]')
legend('Measured', 'Estimate')

% Detailed plot

N2 = 100;
M1_t = Metafile1(1:N2,M_col); % Measured
E1_t = -Metafile1(1:N2,E_col); % Estimate
M2_t = Metafile2(1:N2,M_col);
E2_t = -Metafile2(1:N2,E_col);
M3_t = Metafile3(1:N2,M_col);
E3_t = -Metafile3(1:N2,E_col);
M4_t = Metafile4(1:N2,M_col);
E4_t = -Metafile4(1:N2,E_col);
t = [0:N2-1]'*Ts; % time elapsed

% Overall plot
figure
subplot(2,2,1),plot(t, [M1_t, E1_t], 'LineWidth', 1)
title('Measured vs estimate position plot Q/R=0.001')
grid on
axis tight
xlabel('t [s]')
ylabel('Pos [m]')
legend('Measured', 'Estimate')

subplot(2,2,2),plot(t, [M2_t, E2_t], 'LineWidth', 1)
title('Measured vs estimate position plot Q/R=0.01')
grid on
axis tight
xlabel('t [s]')
ylabel('Pos [m]')
legend('Measured', 'Estimate')

subplot(2,2,3),plot(t, [M3_t, E3_t], 'LineWidth', 1)
title('Measured vs estimate position plot Q/R=0.001 L=-0.002')
grid on
axis tight
xlabel('t [s]')
ylabel('Pos [m]')
legend('Measured', 'Estimate')

subplot(2,2,4),plot(t, [M4_t, E4_t], 'LineWidth', 1)
title('Measured vs estimate position plot Q/R=0.01 L=-0.002')
grid on
axis tight
xlabel('t [s]')
ylabel('Pos [m]')
legend('Measured', 'Estimate')
%display ('Press any key to continue '); pause;