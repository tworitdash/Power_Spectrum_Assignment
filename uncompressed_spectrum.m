% spectrum estimation according to uncompressed signal with noise
clear
load('lowpasssignal.mat');
% load('highpasssignal.mat');
% load('multibandsignal.mat');
% load('bandpasssignal.mat');

%a quick check if the mean of the signal is zero
sum = 0;
for i=1:length(x_sig_and_noise)
   sum = sum + x_sig_and_noise(i); 
end
mean = sum / length(x_sig_and_noise);
display(mean)
x_sig_and_noise = x_sig_and_noise - mean;

% 1. simply compute the periodogram on all samples
Xn = fft(x_sig_and_noise);
Px = Xn.*conj(Xn);
Px = Px/length(Xn);
freq = 0:2*pi/length(Xn):2*pi;
freq = freq(1:length(freq)-1);
freq = freq';
figure()
plot(freq, Px)
set(gca,'XTick',[0 pi/2 pi 3*pi/2 2*pi])
set(gca,'XTickLabel',char('0', '\pi/2', '\pi', '3\pi/2','2\pi'))
xlim([0 2*pi])
ylim([0 60])
xlabel('Radial frequency \omega')
ylabel('P_x(e^j^\omega)')
title('Power Spectrum Estimation Based on the Uncompressed Data with Noise')
hold on

%2. use averaging according to Bartlett's Method
Xk = zeros(length(N_block) ,1);
for l = 1:fix(length(x_sig_and_noise)/N_block)-1
   xl = x_sig_and_noise(l*N_block+1 : (l+1)*N_block);
   Xl = fft(xl);
   Xl = Xl.*conj(Xl);
   Xk = Xk + Xl;
end
Xk = Xk/length(x_sig_and_noise);
freq2 = 0:2*pi/length(Xk):2*pi;
freq2 = freq2(1:length(freq2)-1);
freq2 = freq2';
plot(freq2, Xk);

%3. use minimum variance spectral estimate
K = fix(length(x_sig_and_noise)/N_block);
Rx = zeros(N_block);
for i=0:K-1
   xi = x_sig_and_noise(i*N_block+1:(i+1)*N_block);
   Rx = Rx + xi * xi'; 
end
Rx = Rx/K;
omega = 0:2*pi/N_block:2*pi;
omega = omega(1:length(omega)-1);
omega = omega';
e = NaN(N_block, 1);
Pmv = NaN(length(omega), 1);
for c1=1:length(omega)
    e(1) = 1;
for c2=2:N_block
    e(c2) = exp(1i*(c2-1)*omega(c1));
end
denom = e'/Rx*e;
denom = real(denom);
Pmv(c1) = N_block/denom;
end
plot(omega, Pmv)
legend('Periodogram', 'Bartletts Method', 'Minimum Variance')

% spectrum estimation according to uncompressed signal without noise

%a quick check if the mean of the signal is zero
sum = 0;
for i=1:length(x_sig)
   sum = sum + x_sig(i); 
end
mean = sum / length(x_sig);
display(mean)
x_sig = x_sig - mean;

% 1. simply compute the periodogram on all samples
Xn = fft(x_sig);
Px = Xn.*conj(Xn);
Px = Px/length(Xn);
freq = 0:2*pi/length(Xn):2*pi;
freq = freq(1:length(freq)-1);
freq = freq';
figure()
plot(freq, Px)
set(gca,'XTick',[0 pi/2 pi 3*pi/2 2*pi])
set(gca,'XTickLabel',char('0', '\pi/2', '\pi', '3\pi/2','2\pi'))
xlim([0 2*pi])
ylim([0 60])
xlabel('Radial frequency \omega')
ylabel('P_x(e^j^\omega)')
title('Power Spectrum Estimation Based on the Uncompressed Data without Noise')
hold on

%2. use averaging according to Bartlett's Method
Xk = zeros(length(N_block) ,1);
for l = 1:fix(length(x_sig)/N_block)-1
   xl = x_sig(l*N_block+1 : (l+1)*N_block);
   Xl = fft(xl);
   Xl = Xl.*conj(Xl);
   Xk = Xk + Xl;
end
Xk = Xk/length(x_sig);
freq2 = 0:2*pi/length(Xk):2*pi;
freq2 = freq2(1:length(freq2)-1);
freq2 = freq2';
plot(freq2, Xk);

%3. use minimum variance spectral estimate
K = fix(length(x_sig)/N_block);
Rx = zeros(N_block);
for i=0:K-1
   xi = x_sig(i*N_block+1:(i+1)*N_block);
   Rx = Rx + xi * xi'; 
end
Rx = Rx/K;
omega = 0:2*pi/N_block:2*pi;
omega = omega(1:length(omega)-1);
omega = omega';
e = NaN(N_block, 1);
Pmv = NaN(length(omega), 1);
for c1=1:length(omega)
    e(1) = 1;
for c2=2:N_block
    e(c2) = exp(1i*(c2-1)*omega(c1));
end
denom = e'/Rx*e;
denom = real(denom);
Pmv(c1) = N_block/denom;
end
plot(omega, Pmv)
legend('Periodogram', 'Bartletts Method', 'Minimum Variance')