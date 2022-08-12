clear all
close all
clear global
clc

%% Constants
fs = 100; % sampling frequency (control function - 100Hz)
Ts = 1/fs;

%% Data input and plot
w_col = 8; % supposed to be the same columns for rot.velocity and voltage
v_col = 10; % re-check after each data save
datafile1 = 'SimpleSignal.csv';
datafile2 = 'CompoundSignal.csv';
% -------------------------------------------
Metafile1 = csvread(datafile1); % v_t, w_t, t and nop
Metafile2 = csvread(datafile2); % the same for the compound signal
N = length(Metafile1); % files are supposed to have the same time length ==> check it
v1_t = Metafile1(1:N,v_col); % input voltage
w1_t = Metafile1(1:N,w_col); % rotational velocity
v2_t = Metafile2(1:N,v_col); % input voltage
w2_t = Metafile2(1:N,w_col); % rotational velocity

t = [0:N-1]'*Ts; % time elapsed
nop = 5; % (adjust if needed according to data sample)
ppp = N/nop;
figure
subplot(2,1,1),plot(t, w1_t, t, w2_t, 'LineWidth', 1)
title('Rotational velocity plot')
grid on, axis tight
xlabel('t [s]')
ylabel('omega [rad/s]')
subplot(2,1,2), plot(t, v1_t, t, v2_t, 'LineWidth', 1)
title('Voltage profile')
grid on, axis tight
xlabel('t [s]')
ylabel('Voltage [V]')
w_matrix1 = reshape(w1_t,ppp,nop);
w_matrix1 = w_matrix1 - ones(ppp,1)*w_matrix1(1,:);
dw_matrix1 = w_matrix1 - [w_matrix1(:,2:end), w_matrix1(:,1)];
v_matrix1 = reshape(v1_t,ppp,nop);
dv_matrix1 = v_matrix1 - [v_matrix1(:,2:end), v_matrix1(:,1)];
w_matrix2 = reshape(w2_t,ppp,nop);
w_matrix2 = w_matrix2 - ones(ppp,1)*w_matrix2(1,:);
dw_matrix2 = w_matrix2 - [w_matrix2(:,2:end), w_matrix2(:,1)];
v_matrix2 = reshape(v2_t,ppp,nop);
dv_matrix2 = v_matrix2 - [v_matrix2(:,2:end), v_matrix2(:,1)];
figure,subplot(2,2,1),plot(t(1:ppp), w_matrix1, t(1:ppp), w_matrix2, 'LineWidth', 1)
grid on, axis tight
xlabel('t  [s]')
ylabel('omega  [rad/s]')
subplot(2,2,3),plot(t(1:ppp), v_matrix1, t(1:ppp), v_matrix2, 'LineWidth', 1)
grid on, axis tight
xlabel('t  [s]')
ylabel('Voltage  [V]')
subplot(2,2,2),plot(t(1:ppp), dw_matrix1, t(1:ppp), dw_matrix2, 'LineWidth', 1)
grid on, axis tight
xlabel('t  [s]')
ylabel('\Delta omega  [rad/s]')
subplot(2,2,4)
plot(t(1:ppp), dv_matrix1, t(1:ppp), dv_matrix2, 'LineWidth', 1)
grid on, axis tight
xlabel('t  [s]')
ylabel('\Delta V  [V]')
%display ('Press any key to continue '); pause;

%% Calculating and plotting the empirical transfer-functions estimate 
N = numel(v1_t);
f = [0:N-1]'*(fs/N);
v1_f = fft(v1_t);
w1_f = fft(w1_t);
v2_f = fft(v2_t);
w2_f = fft(w2_t);
indices = (nop+1):nop:(numel(f)/2);
f = f(indices);
v1_f = v1_f(indices);
w1_f = w1_f(indices);
v2_f = v2_f(indices);
w2_f = w2_f(indices);
FRF1 = w1_f./v1_f; % Empirical FRF - original signal, further FRF1_1, FRF1_2 = its LLS estimations
FRF2 = w2_f./v2_f; % Empirical FRF - compound signal, further FRF2_1, FRF2_2 = its LLS estimations
figure,subplot(2,1,1),semilogx(f, 20*log10(abs(FRF1)), f, 20*log10(abs(FRF2)), 'LineWidth', 1)
title('Bode plot of the empirical transfer functions (magnitude)')
grid on, axis tight
xlim([f(1) f(end)])
xlabel('f [Hz]')
ylabel('|FRF| [dB]')
legend({'Simple signal', 'Compound signal'}, 'Location', 'northwest')
subplot(2,1,2),semilogx(f, 180/pi*unwrap(angle(FRF1)), f, 180/pi*unwrap(angle(FRF2)), 'LineWidth', 1)
title('Bode plot of the empirical transfer functions (phase)')
grid on, axis tight
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('\phi(FRF) [^\circ]')
legend({'Simple signal', 'Compound signal'}, 'Location', 'southwest')
%display ('Press any key to continue '); pause;

%% LLS without data filtering
% -------------------------------------------
B1 = w1_t(2:end);
A1 = [-w1_t(1:end-1), v1_t(2:end)];
theta1 = A1\B1;
B1_1 = [theta1(2), 0];
A1_1 = [1, theta1(1)];
sys_d1_1 = tf(B1_1, A1_1, Ts);
FRF1_1 = squeeze(freqresp(sys_d1_1,2*pi*f));

B2 = w2_t(2:end);
A2 = [-w2_t(1:end-1), v2_t(2:end)];
theta2 = A2\B2;
B2_1 = [theta2(2), 0];
A2_1 = [1, theta2(1)];
sys_d2_1 = tf(B2_1, A2_1, Ts);
FRF2_1 = squeeze(freqresp(sys_d2_1,2*pi*f));

figure
subplot(2,1,1)
semilogx(f, 20*log10(abs(FRF1)), f, 20*log10(abs(FRF2)), f, 20*log10(abs(FRF1_1)), f, 20*log10(abs(FRF2_1)), 'LineWidth', 1);
title('Bode plot of empirical/LLS TF (magnitude)')
grid on
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('|FRF|  [dB]')
legend({'Empirical Simple', 'Empirical Comp', 'Unfiltered Simple', 'Unfiltered Comp'}, 'Location', 'northwest')
subplot(2,1,2)
semilogx(f, 180/pi*unwrap(angle(FRF1)), f, 180/pi*unwrap(angle(FRF2)), f, 180/pi*unwrap(angle(FRF1_1)), f, 180/pi*unwrap(angle(FRF2_1)), 'LineWidth', 1);
title('Bode plot of empirical/LLS TF (phase)')
grid on
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('\phi(FRF)  [^\circ]')
legend({'Empirical Simple', 'Empirical Comp', 'Unfiltered Simple', 'Unfiltered Comp'}, 'Location', 'southwest')
w1_1 = lsim(sys_d1_1,v1_t,t);
w2_1 = lsim(sys_d2_1,v2_t,t);

figure
subplot(211),plot(t, [w1_t w1_1], t, [w2_t w2_1]);
legend('w simp emp', 'w simp recon', 'w comp emp', 'w comp recon')
xlabel('time [s]')
ylabel('rotational velocity [rad/s]')
axis tight
subplot(212),plot(t, w1_t - w1_1, t, w2_t - w2_1)
legend('w simp difference', 'w comp difference')
xlabel('time [s]')
ylabel('rotational velocity [rad/s]')
axis tight
figure,pzmap(sys_d1_1, sys_d2_1)
%display ('Press any key to continue '); pause;

%% 3.3.2 LLS with low-pass filter applied to the input and output data
% -------------------------------------------------------------------
[B_filt,A_filt] = butter(6, (2/fs)*10);
w1_filt = filter(B_filt, A_filt, w1_t); 
v1_filt = filter(B_filt, A_filt, v1_t);
w2_filt = filter(B_filt, A_filt, w2_t); 
v2_filt = filter(B_filt, A_filt, v2_t);

B1_1 = w1_filt(2:end);
A1_1 = [-w1_filt(1:end-1), v1_filt(2:end)];
theta1 = A1_1\B1_1;
B2_1 = [theta1(2), 0];
A2_1 = [1, theta1(1)];
sys_d1_2 = tf(B2_1, A2_1, Ts);
FRF1_2 = squeeze(freqresp(sys_d1_2,2*pi*f));

B1_2 = w2_filt(2:end);
A1_2 = [-w2_filt(1:end-1), v2_filt(2:end)];
theta2 = A1_2\B1_2;
B2_2 = [theta2(2), 0];
A2_2 = [1, theta2(1)];
sys_d2_2 = tf(B2_2, A2_2, Ts);
FRF2_2 = squeeze(freqresp(sys_d2_2,2*pi*f));

figure
subplot(2,1,1)
semilogx(f, 20*log10(abs(FRF1)), f, 20*log10(abs(FRF2)), f, 20*log10(abs(FRF1_2)), f, 20*log10(abs(FRF2_2)), 'LineWidth', 1)
title('Bode plot of empirical/LLS TF filt (magnitude)')
grid on,axis tight
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('|FRF|  [dB]')
legend({'Empirical Simple', 'Empirical Comp', 'Filtered Simple', 'Filtered Comp'}, 'Location', 'northwest')
subplot(2,1,2)
semilogx(f, 180/pi*unwrap(angle(FRF1)), f, 180/pi*unwrap(angle(FRF2)), f, 180/pi*unwrap(angle(FRF1_2)), f, 180/pi*unwrap(angle(FRF2_2)), 'LineWidth', 1)
title('Bode plot of empirical/LLS TF filt (phase)')
grid on,axis tight
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('\phi(FRF)  [^\circ]')
legend({'Empirical Simple', 'Empirical Comp', 'Filtered Simple', 'Filtered Comp'}, 'Location', 'southwest')

w1_2 = lsim(sys_d1_2, v1_t, t);
w2_2 = lsim(sys_d2_2, v2_t, t);
figure,subplot(211),plot(t, [w1_t w1_2], t , [w2_t w2_2]);
legend('w simp emp', 'w simp recon', 'w comp emp', 'w comp recon')
xlabel('time [s]')
axis tight
ylabel('rotational velocity [rad/s]')
subplot(212),plot(t, w1_t-w1_2, t, w2_t-w2_2)
legend('w simp difference', 'w comp difference')
xlabel('time [s]')
ylabel('rotational velocity [rad/s]')
axis tight
figure,pzmap(sys_d1_2, sys_d2_2)
%display ('Press any key to continue '); pause;

%% Bode plots (alle functies):
figure
subplot(2,1,1)
semilogx(f, 20*log10(abs(FRF1)), f, 20*log10(abs(FRF1_1)), f, 20*log10(abs(FRF1_2)), f, 20*log10(abs(FRF2)), f, 20*log10(abs(FRF2_1)), f, 20*log10(abs(FRF2_2)), 'LineWidth', 1)
title('Bode plot of transfer functions (magnitude)')
grid on, axis tight
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('|FRF|  [dB]')
legend({'Empirical Simple', 'Unfiltered Simple', 'Filtered Simple', 'Empirical Comp', 'Unfiltered Comp', 'Filtered Comp'}, 'Location', 'northwest')
subplot(2,1,2)
semilogx(f, 180/pi*unwrap(angle(FRF1)), f, 180/pi*unwrap(angle(FRF1_1)), f, 180/pi*unwrap(angle(FRF1_2)), f, 180/pi*unwrap(angle(FRF2)), f, 180/pi*unwrap(angle(FRF2_1)), f, 180/pi*unwrap(angle(FRF2_2)), 'LineWidth', 1)
title('Bode plot of transfer functions (phase)')
grid on, axis tight
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('\phi(FRF)  [^\circ]')
legend({'Empirical Simple', 'Unfiltered Simple', 'Filtered Simple', 'Empirical Comp', 'Unfiltered Comp', 'Filtered Comp'}, 'Location', 'southwest')

%% Bode plots (enkel estimaties):
figure
subplot(2,1,1)
semilogx(f, 20*log10(abs(FRF1_1)), f, 20*log10(abs(FRF1_2)), f, 20*log10(abs(FRF2_1)), f, 20*log10(abs(FRF2_2)), 'LineWidth', 1)
title('Bode plot of estimated TF (magnitude)')
grid on, axis tight
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('|FRF|  [dB]')
legend({'Unfiltered Simple', 'Filtered Simple', 'Unfiltered Comp', 'Filtered Comp'}, 'Location', 'southwest')
subplot(2,1,2)
semilogx(f, 180/pi*unwrap(angle(FRF1_1)), f, 180/pi*unwrap(angle(FRF1_2)), f, 180/pi*unwrap(angle(FRF2_1)), f, 180/pi*unwrap(angle(FRF2_2)), 'LineWidth', 1)
title('Bode plot of estimated TF (phase)')
grid on, axis tight
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('\phi(FRF)  [^\circ]')
legend({'Unfiltered Simple', 'Filtered Simple', 'Unfiltered Comp', 'Filtered Comp'}, 'Location', 'southwest')
 
%% Bode plots (enkel ongefilterde)
figure
subplot(2,1,1)
semilogx(f, 20*log10(abs(FRF1_1)), f, 20*log10(abs(FRF2_1)), 'LineWidth', 1)
title('Bode plot of LLS estimations (magnitude)')
grid on, axis tight
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('|FRF|  [dB]')
legend({'Unfiltered Simple', 'Unfiltered Comp'}, 'Location', 'southwest')
subplot(2,1,2)
semilogx(f, 180/pi*unwrap(angle(FRF1_1)), f, 180/pi*unwrap(angle(FRF2_1)), 'LineWidth', 1)
title('Bode plot of LLS estimations (phase)')
grid on, axis tight
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('\phi(FRF)  [^\circ]')
legend({'Unfiltered Simple', 'Unfiltered Comp'}, 'Location', 'southwest')

%% Bode plots (enkel gefilterde)
figure
subplot(2,1,1)
semilogx(f, 20*log10(abs(FRF1_2)), f, 20*log10(abs(FRF2_2)), 'LineWidth', 1)
title('Bode plot of LLS estimations filtered (magnitude)')
grid on, axis tight
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('|FRF|  [dB]')
legend({'Filtered Simple', 'Filtered Comp'}, 'Location', 'southwest')
subplot(2,1,2)
semilogx(f, 180/pi*unwrap(angle(FRF1_2)), f, 180/pi*unwrap(angle(FRF2_2)), 'LineWidth', 1)
title('Bode plot of LLS estimations filtered (phase)')
grid on, axis tight
xlim([f(1) f(end)])
xlabel('f  [Hz]')
ylabel('\phi(FRF)  [^\circ]')
legend({'Filtered Simple', 'Filtered Comp'}, 'Location', 'southwest')