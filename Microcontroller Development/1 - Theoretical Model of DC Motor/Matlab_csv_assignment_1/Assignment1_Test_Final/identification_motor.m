clear all
close all
clear global
clc

%% Constants
fs = 100; % sampling frequency (control function - 100Hz)
Ts = 1/fs;

%% Data input and plot
w_col = 7;
v_col = 9;
datafile = 'leftWheelFinal.csv';
% -------------------------------------------
Metafile = csvread(datafile); % v_t, w_t, t and nop
N = length(Metafile);
v_t = Metafile(1:N,v_col); % input voltage
w_t = Metafile(1:N,w_col); % rotational velocity
t = [0:N-1]'*Ts; % time elapsed
nop = 5; % (adjust if needed according to data sample)
ppp = N/nop;
figure,subplot(2,1,1),plot(t, w_t, 'LineWidth', 1)
title('Rotational velocity plot')
grid on
axis tight
xlabel('t [s]')
ylabel('omega [rad/s]')
subplot(2,1,2), plot(t, v_t, 'LineWidth', 1)
title('Voltage profile')
grid on
axis tight
xlabel('t [s]')
ylabel('Voltage [V]')
w_matrix = reshape(w_t,ppp,nop);
w_matrix = w_matrix - ones(ppp,1)*w_matrix(1,:);
dw_matrix = w_matrix - [w_matrix(:,2:end), w_matrix(:,1)];
v_matrix = reshape(v_t,ppp,nop);
dv_matrix = v_matrix - [v_matrix(:,2:end), v_matrix(:,1)];
figure,subplot(2,2,1),plot(t(1:ppp), w_matrix, 'LineWidth', 1)
grid on
axis tight
xlabel('t  [s]')
ylabel('omega  [rad/s]')
subplot(2,2,3),plot(t(1:ppp), v_matrix, 'LineWidth', 1)
grid on
axis tight
xlabel('t  [s]')
ylabel('Voltage  [V]')
subplot(2,2,2),plot(t(1:ppp), dw_matrix, 'LineWidth', 1)
grid on
axis tight
xlabel('t  [s]')
ylabel('\Delta omega  [rad/s]')
subplot(2,2,4)
plot(t(1:ppp), dv_matrix, 'LineWidth', 1)
grid on
axis tight
xlabel('t  [s]')
ylabel('\Delta V  [V]')
%display ('Press any key to continue '); pause;

%% Calculating and plotting the empirical transfer-function estimate 
N = numel(v_t);
f = [0:N-1]'*(fs/N);
v_f = fft(v_t);
w_f = fft(w_t);
indices = (nop+1):nop:(numel(f)/2);
f = f(indices);
v_f = v_f(indices);
w_f = w_f(indices);
FRF = w_f./v_f; % Empirical FRF
axis tight
figure
subplot(2,1,1),semilogx(f, 20*log10(abs(FRF)), 'LineWidth', 1)
title('Bode plot of the empirical transfer function (magnitude)')
grid on
xlabel('f [Hz]')
xlim([f(1) f(end)])
ylabel('|FRF| [dB]')
subplot(2,1,2),semilogx(f, 180/pi*unwrap(angle(FRF)), 'LineWidth', 1)
title('Bode plot of the empirical transfer function (phase)')
grid on
axis tight
xlabel('f  [Hz]')
ylabel('\phi(FRF) [^\circ]')
xlim([f(1) f(end)])
%display ('Press any key to continue '); pause;

%% LLS without data filtering
% -------------------------------------------
B = w_t(2:end);
A = [-w_t(1:end-1), v_t(2:end)];
theta = A\B;
B1 = [theta(2), 0];
A1 = [1, theta(1)];
sys_d1 = tf(B1, A1, Ts);
FRF1 = squeeze(freqresp(sys_d1,2*pi*f));
figure
subplot(2,1,1),semilogx(f, 20*log10(abs(FRF)), f, 20*log10(abs(FRF1)))
title('Bode plot of transfer functions (magnitude)')
grid on
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('|FRF|  [dB]')
legend('Empirical', 'Estimated')
subplot(2,1,2),semilogx(f, 180/pi*unwrap(angle(FRF)), f, 180/pi*unwrap(angle(FRF1)))
title('Bode plot of transfer functions (phase)')
grid on
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('\phi(FRF)  [^\circ]')
legend('Empirical', 'Estimated')
w1 = lsim(sys_d1,v_t,t);
figure,subplot(211),plot(t,[w_t w1]);
legend('w_t','w_1')
xlabel('time [s]')
ylabel('rotational velocity [rad/s]')
axis tight
subplot(212),plot(t,w_t - w1)
legend('w_t-w_1')
xlabel('time [s]')
ylabel('rotational velocity [rad/s]')
axis tight
figure,pzmap(sys_d1)
%display ('Press any key to continue '); pause;

%% 3.3.2 LLS with low-pass filter applied to the input and output data
% -------------------------------------------------------------------
[B_filt,A_filt] = butter(6, (2/fs)*10);
w_filt = filter(B_filt, A_filt, w_t); 
v_filt = filter(B_filt, A_filt, v_t);
B = w_filt(2:end);
A = [-w_filt(1:end-1), v_filt(2:end)];
theta = A\B;
B2 = [theta(2), 0];
A2 = [1, theta(1)];
sys_d2 = tf(B2, A2, Ts);
FRF2 = squeeze(freqresp(sys_d2,2*pi*f));
figure
subplot(2,1,1),semilogx(f, 20*log10(abs(FRF)), f, 20*log10(abs(FRF2)))
title('Bode plot of transfer functions (magnitude)')
grid on
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('|FRF|  [dB]')
legend('Empirical', 'Filtered')
axis tight
subplot(2,1,2),semilogx(f, 180/pi*unwrap(angle(FRF)), f, 180/pi*unwrap(angle(FRF2)))
title('Bode plot of transfer functions (phase)')
grid on
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('\phi(FRF)  [^\circ]')
legend('Empirical', 'Filtered')
w2 = lsim(sys_d2,v_t,t);
figure,subplot(211),plot(t,[w_t w2]);
legend('w_s','w_f')
xlabel('time [s]')
axis tight
ylabel('rotational velocity [rad/s]')
subplot(212),plot(t,w_t-w2)
legend('w_s-w_f')
xlabel('time [s]')
ylabel('rotational velocity [rad/s]')
axis tight
figure,pzmap(sys_d2)
%display ('Press any key to continue '); pause;

%% Bode plots:
figure
subplot(2,1,1)
semilogx(f, 20*log10(abs(FRF)), f, 20*log10(abs(FRF1)), f, 20*log10(abs(FRF2)))
title('Bode plot of transfer functions (magnitude)')
grid on
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('|FRF|  [dB]')
legend('Empirical', 'Unfiltered', 'Filtered')
axis tight
subplot(2,1,2)
semilogx(f, 180/pi*unwrap(angle(FRF)), f, 180/pi*unwrap(angle(FRF1)), f, 180/pi*unwrap(angle(FRF2)))
title('Bode plot of transfer functions (phase)')
grid on
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('\phi(FRF)  [^\circ]')
legend('Empirical', 'Unfiltered', 'Filtered')

% Die twee samen
figure
subplot(2,1,1)
semilogx(f, 20*log10(abs(FRF1)), f, 20*log10(abs(FRF2)))
title('Bode plot of transfer functions (magnitude)')
grid on
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('|FRF|  [dB]')
legend('Unfiltered', 'Filtered')
axis tight
subplot(2,1,2)
semilogx(f, 180/pi*unwrap(angle(FRF1)), f, 180/pi*unwrap(angle(FRF2)))
title('Bode plot of transfer functions (phase)')
grid on
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('\phi(FRF)  [^\circ]')
legend('Unfiltered', 'Filtered')