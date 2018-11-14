
%Computation of Phi (The matrix of N_block and M_ruler) to compress the
%exisiting data.
load('lowpasssignal.mat')
mn = dsp.Mean;
y_sig_mean = mn(y_sig);
y_sig_updated = y_sig - y_sig_mean;

Phi = zeros(M_ruler, N_block);
for i = 1:M_ruler-1
    Phi(i, ruler_index(i)) = 1;
end;

%computation of Ry - Autocorrelation of y
