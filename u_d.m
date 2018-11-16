%spectrum of uncompressed signal x
load('lowpasssignal.mat')
%check for mean
mn = dsp.Mean;
x_sig_with_noise_mean = mn(x_sig_and_noise);

%new signal definition with mean

x_sig_with_noise_updated = x_sig_and_noise - x_sig_with_noise_mean;

%Periodogram on x to compute PSD
N = length(x_sig_with_noise_updated);
