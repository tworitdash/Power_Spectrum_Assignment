%spectrum_uncompressed.m

%different mat files can be used by just removing the comment below. Now, it is active for lowpass.mat

%spectrum of uncompressed signal x

load('lowpasssignal.mat');
%load('multibandsignal.mat');
%load('highpasssignal.mat');
%load('bandpasssignal.mat');

%check for mean

mn = dsp.Mean;
x_sig_with_noise_mean = mn(x_sig_and_noise);

%new signal definition with mean

x_sig_with_noise_updated = x_sig_and_noise - x_sig_with_noise_mean;

%Periodogram on x to compute PSD

N = length(x_sig_with_noise_updated);

X_per = periodo(x_sig_with_noise_updated, N);

fs = 2*pi; 

%sampling of 2*pi becuase signal is complex valued
    
freq_0 = 0:fs/length(x_sig_with_noise_updated):fs - fs/N;
figure(1)

plot(freq_0/pi, abs(X_per));
grid on
title('PSD with Periodogram using FFT and Bartlett Band Pass')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Power/Frequency')

%Bartlett's method

X_per_bartlett = Pxl(x_sig_with_noise_updated, N_block);

freq_1 = 0:fs/length(X_per_bartlett):fs;
freq_1_updated = freq_1(1:length(freq_1)-1);

hold on
plot(freq_1_updated/pi, (X_per_bartlett), 'LineWidth',3);
grid on
%title('PSD with Periodogram using Bartlett Avg')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Power/Frequency')



legend('FFT', 'Bartlett')

%function periodo for general periodogram

function perx = periodo(x, N)
  
        x_fft = fft(x);
        perx = (1/ N) * (x_fft) .* conj(x_fft);
    
end


%function Pxl for Bartlett

function per_x = Pxl(x, N_block)

    X = zeros(length(N_block) ,1);

    for i = 1:fix(length(x)/N_block)-1 %finding it out for different sections
        
       x_i = x(i * N_block + 1 : (i + 1) * N_block);

       X_i = fft(x_i);

       X_i = abs(X_i) .^ 2;

       X = X + X_i;
    end
    
    per_x = X/length(x);
end
