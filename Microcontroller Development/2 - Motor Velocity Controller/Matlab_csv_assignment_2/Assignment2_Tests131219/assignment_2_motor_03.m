clear all
close all
clear global
clc

fs = 100;
Ts = 1/fs;
%% Import of data from the experiment

% Registered data:
% 2-3: left/right rotational velocity
% 4  : reference (desired) velocity
% 5-6: left/right velocity error
% 7-8: left/right engine voltage
w_col = 2;
wr_col = 4;
e_col = 5;
v_col = 7;
datafile = 'Control03Wheel65.csv';
% -------------------------------------------
Metafile = csvread(datafile); % v_t, w_t, t and nop
N = length(Metafile); % files are supposed to have the same time length ==> check it
w_exp = Metafile(1:N,w_col);
wr_exp = Metafile(1:N,wr_col);
e_exp = Metafile(1:N,e_col);
v_exp = Metafile(1:N,v_col);

t = [0:N-1]'*Ts; % time elapsed

nop = 5; % (adjust if needed according to data sample)
ppp = N/nop;

% No plot of experimental function (add later if necessary)
figure
subplot(2,1,1),plot(t, w_exp, 'LineWidth', 1)
title('Rotational velocity profile')
grid on
axis tight
xlabel('t [s]')
ylabel('omega [rad/s]')
subplot(2,1,2), plot(t, v_exp, 'LineWidth', 1)
title('Voltage profile')
grid on
axis tight
xlabel('t [s]')
ylabel('Voltage [V]')

%% Continous function of the motor:

% Discrete function based on motor experiment

b1_d = 0.3970;
a1_d = 1;
a0_d = -0.8109;

% Manual conversion of motor function to continuous with Backward Euler
b0_c = -b1_d/(a0_d*Ts);
a1_c = 1;
a0_c = -(a1_d + a0_d)/(a0_d*Ts);

num_G = b0_c;
den_G = [a1_c a0_c];

sys_G = tf(num_G, den_G);

% 1. Motor TF Bode-diagram
w = logspace(-0.3,1.7,100);		%100 frequentiepoints, logarithmically spaced between 10^0 and 10^2 rad/s
[mag_G, phase_G] = bode(num_G,den_G,w);

% Bode diagram of the continuous transfer function
figure
subplot(211),semilogx(w, 20*log10(mag_G))
title('Bodeplot (magnitude) of the initial TF')
grid on, axis tight
xlabel('w [rad/s]')
ylabel('|G(j\omega)| [dB]')
subplot(212),semilogx(w,phase_G)
title('Bodeplot (phase) of the initial TF')
grid on, axis tight
xlabel('w [rad/s]')
ylabel('\phi(G(j\omega)) [^o]')
 
%% Controller design

% 2.1. Determine the new cross-over pulsation wco:
Dphi_PI = 15; % Desired phase lag (degrees)
wco = 10; % Crossover frequency

% 2.2. Determine Ti [s], such that the phase lag of the PI controller at wco equals Dphi_PI
%Ti = 1/(wco * tan(Dphi_PI*pi/180));      % Reduced due to our needs (period of 1s)
Ti = 0.3732; % adapt when necessary

 % PI controller parameters
num_PI = [Ti 1];
den_PI = [Ti 0];

% 2.3. Calculate the gain such that the amplitude at wco equals 1
num_open = conv(num_G, num_PI);
den_open = conv(den_G, den_PI);
[mag_open, phase_open] = bode(num_open, den_open, w);
%gain = 1/interp1(w,mag_loop,wco); %Changed ==> adapt according to experim
gain = 0.4713; % adapt when necessary

% 2.4. Bodeplot of the compensated system (OPEN LOOP)
num_PI = num_PI*gain;
sys_PI = tf(num_PI, den_PI);
num_open = num_open*gain;
sys_open = tf(num_open, den_open);
[mag_open,phase_open] = bode(num_open, den_open, w);

% Open-loop
figure
subplot(211),semilogx(w, 20*log10(mag_open))
title('Bodeplot (magnitude) of the PI-compensated transfer function (open loop)')
grid on, axis tight
xlabel('w [rad/s]')
ylabel('|G(j\omega)| [dB]')
subplot(212),semilogx(w, phase_open)
title('Bodeplot (phase) of the PI-compensated transfer function (open loop)')
grid on, axis tight
xlabel('w [rad/s]')
ylabel('\phi(G(j\omega)) [^o]')

% Closed-loop response transfer function
sys_closed = feedback(sys_open, 1);
[mag_closed,phase_closed] = bode(sys_closed, w);
mag_closed = squeeze(mag_closed);
phase_closed = squeeze(phase_closed);

% Closed-loop
figure
subplot(211),semilogx(w, 20*log10(mag_closed))
title('Bodeplot (magnitude) of the PI-compensated transfer function (closed loop)')
grid on, axis tight
xlabel('w [rad/s]')
ylabel('|G(j\omega)| [dB]')
subplot(212),semilogx(w, phase_closed)
title('Bodeplot (phase) of the PI-compensated transfer function (closed loop)')
grid on, axis tight
xlabel('w [rad/s]')
ylabel('\phi(G(j\omega)) [^o]')

 % Uncompensated vs open vs closed (all continuous)
figure
subplot(211),semilogx(w, 20*log10([mag_G, mag_open, mag_closed]), 'Linewidth', 1)
title('Bodeplot (magnitude)')
grid on, axis tight
xlabel('w [rad/s]')
ylabel('|G(j\omega)| [dB]')
legend('Uncompensated', 'Open loop', 'Closed loop')
subplot(212), semilogx(w, [phase_G, phase_open, phase_closed], 'Linewidth', 1)
title('Bodeplot (phase)')
grid on, axis tight
xlabel('w [rad/s]')
ylabel('\phi(G(j\omega)) [^o]')
legend('Uncompensated', 'Open loop', 'Closed loop')

%% Verification of the chain
% Application of our continous transfer function to velocity data

% % Check 1: Voltage to error relation (PI-regelaar)
% v_th = lsim(sys_PI,e_exp,t);
% figure
% subplot(211),plot(t,[v_exp v_th]);
% title('Check of the PI-controller')
% grid on, axis tight
% xlabel('time [s]')
% ylabel('Voltage [V]')
% legend('Voltage exper','Voltage based on TF')
% subplot(212),plot(t,v_exp - v_th)
% grid on, axis tight
% legend('Difference exper - based on TF')
% xlabel('time [s]')
% ylabel('Voltage [V]')
% 
% % Check 2: Velocity to voltage relation ("Original" transfer function)
% w_th = lsim(sys_G,v_exp,t);
% figure
% subplot(211),plot(t,[w_exp w_th]);
% title('Check of the transfer function')
% grid on, axis tight
% xlabel('time [s]')
% ylabel('w [rad/s]')
% legend('Rot.velocity exper','Rot.velocity based on TF')
% subplot(212),plot(t,w_exp - w_th)
% title('Check of the transfer function')
% grid on, axis tight
% xlabel('time [s]')
% ylabel('w [rad/s]')
% legend('Rot.velocity exper - based on TF')
% 
% % Check 3: Velocity to error relation (open loop)
% w_th = lsim(sys_open,e_exp,t);
% figure
% subplot(211),plot(t,[w_exp w_th]);
% title('Check of the open loop')
% grid on, axis tight
% xlabel('time [s]')
% ylabel('w [rad/s]')
% legend('Rot.velocity exper','Rot.velocity based on open loop')
% subplot(212),plot(t,w_exp - w_th)
% title('Check of the open loop')
% grid on, axis tight
% xlabel('time [s]')
% ylabel('w [rad/s]')
% legend('Rot.velocity exper - based on open loop')


%% Task 2a - validate the controller experimentally

% Check 1: Desired vs measured vs simulated velocity
w_th = lsim(sys_closed, w_exp, t);
figure
subplot(211),plot(t,[wr_exp w_exp w_th], 'Linewidth', 1);
title('Check 1: Desired vs measured vs simulated velocity')
grid on, axis tight
xlabel('time [s]')
ylabel('w [rad/s]')
legend('Desired velocity', 'Measured CL response', 'Simulated CL response')
subplot(212),plot(t,w_exp - w_th, 'Linewidth', 1)
title('Difference measured - simulated velocity')
grid on, axis tight
xlabel('time [s]')
ylabel('w [rad/s]')
legend('Meas-sim CL response')

% Check 2: Measured tracking error vs simulated tracking error
err_closed = feedback(1, sys_open);
e_th = lsim(err_closed, wr_exp, t);
figure
subplot(211),plot(t,[e_exp e_th], 'Linewidth', 1);
title('Check 2: Measured vs simulated tracking error')
grid on, axis tight
xlabel('time [s]')
ylabel('w [rad/s]')
legend('Measured track error', 'Simulated track error')
subplot(212),plot(t,e_exp - e_th, 'Linewidth', 1)
title('Difference measured - simulated tracking error')
grid on, axis tight
xlabel('time [s]')
ylabel('w [rad/s]')
legend('Meas-sim tracking error')

% Check 3: Measured control signal vs simulated control signal
 volt_closed = feedback(sys_PI, sys_G);
 v_th = lsim(volt_closed, wr_exp, t);
%v_th = lsim(sys_PI, e_exp, t); % eigenlijk the same
figure
subplot(211),plot(t,[v_exp v_th], 'Linewidth', 1);
title('Check 3: Measured vs simulated control signal')
grid on, axis tight
xlabel('time [s]')
ylabel('Control signal [V]')
legend('Measured control signal', 'Simulated control signal')
subplot(212),plot(t,v_exp - v_th, 'Linewidth', 1)
title('Difference measured - simulated control signal')
grid on, axis tight
xlabel('time [s]')
ylabel('Control signal [V]')
legend('Meas-sim control signal')