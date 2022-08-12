clear all
close all
clear global
clc

%% Constants
fs = 100; % sampling frequency (control function - 100Hz)
Ts = 1/fs;

%% Data input and plot
P_col = 11;
L_col = 12;
datafile1 = 'PLtestQR0001.csv'; % csvread werkt niet echt met matrices
datafile2 = 'PLtestQR001.csv';
datafile3 = 'PLtestQR01.csv';
datafile4 = 'PLtestQR1.csv';
Metafile1 = csvread(datafile1); % P_t and L_t
Metafile2 = csvread(datafile2); % P_t and L_t
Metafile3 = csvread(datafile3); % P_t and L_t
Metafile4 = csvread(datafile4); % P_t and L_t
% -------------------------------------------
N = 1000; % adapted to first 10s of the signal
P1_t = Metafile1(1:N,P_col); % Pkk
L1_t = Metafile1(1:N,L_col); % Lkk
P2_t = Metafile2(1:N,P_col); % Pkk
L2_t = Metafile2(1:N,L_col); % Lkk
P3_t = Metafile3(1:N,P_col); % Pkk
L3_t = Metafile3(1:N,L_col); % Lkk
P4_t = Metafile4(1:N,P_col); % Pkk
L4_t = Metafile4(1:N,L_col); % Lkk
t = [0:N-1]'*Ts; % time elapsed
%nop = 5; % (adjust if needed according to data sample)
%ppp = N/nop;
figure,subplot(2,1,1),plot(t, [P1_t, P2_t, P3_t, P4_t], 'LineWidth', 1)
title('State estimate covariance (P) plot')
grid on
axis tight
xlabel('t [s]')
ylabel('P [m^2]')
legend('Q/R=0.001', 'Q/R=0.01', 'Q/R=0.1', 'Q/R=1')
subplot(2,1,2), plot(t, [L1_t, L2_t, L3_t, L4_t], 'LineWidth', 1)
title('Kalman gain (L) plot')
grid on
axis tight
xlabel('t [s]')
ylabel('L [-]')
legend('Q/R=0.001', 'Q/R=0.01', 'Q/R=0.1', 'Q/R=1')

% Detailed plot

N2 = 50;
P1_t = Metafile1(1:N2,P_col); % Measured
L1_t = Metafile1(1:N2,L_col); % Estimate
P2_t = Metafile2(1:N2,P_col);
L2_t = Metafile2(1:N2,L_col);
P3_t = Metafile3(1:N2,P_col);
L3_t = Metafile3(1:N2,L_col);
P4_t = Metafile4(1:N2,P_col);
L4_t = Metafile4(1:N2,L_col);
t = [0:N2-1]'*Ts; % time elapsed

% Overall plot
figure
subplot(2,1,1),plot(t, [P1_t, P2_t, P3_t, P4_t], 'LineWidth', 1)
title('State estimate covariance (P) plot')
grid on
axis tight
xlabel('t [s]')
ylabel('P [m^2]')
legend('Q/R=0.001', 'Q/R=0.01', 'Q/R=0.1', 'Q/R=1')

subplot(2,1,2),plot(t, [L1_t, L2_t, L3_t, L4_t], 'LineWidth', 1)
title('Kalman gain (L) plot')
grid on
axis tight
xlabel('t [s]')
ylabel('L [-]')
legend('Q/R=0.001', 'Q/R=0.01', 'Q/R=0.1', 'Q/R=1')

%display ('Press any key to continue '); pause;