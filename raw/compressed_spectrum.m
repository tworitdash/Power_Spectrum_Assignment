%compression according to the master thesis Chapter 2
clear
load('lowpasssignal.mat');
% load('highpasssignal.mat');
% load('multibandsignal.mat');
% load('bandpasssignal.mat');

%a quick check if the mean of the signal is zero
sum = 0;
for i = 1:length(y_sig)
   sum = sum + y_sig(i); 
end
mean = sum/length(y_sig);
display(mean)
y_sig = y_sig - mean;

%compute sampling matrix C
C = zeros(M_ruler, N_block);
for i=1:M_ruler
   C(i,ruler_index(i)) = 1; 
end

%compute the autocorrelation matrix of y_sig according to the formulas
%given in the master thesis and the paper spawk.pdf
K = length(x_sig_and_noise)/N_block;

VecRy = zeros(M_ruler^2,1);
for i=1:M_ruler
   ryi0 = zeros(M_ruler, 1);
   for j = 1:M_ruler
      ryi0(j) = ry(j-1, i-1, y_sig, M_ruler);
   end
   VecRy((i-1)*M_ruler+1:i*M_ruler) = ryi0;
end

%compute T according to the formula given in the master thesis
T = zeros(N_block^2, 2*N_block-1);
Iaux = eye(2*N_block-1);
for i = 1:N_block^2
   T(i,:) = Iaux(mod(i-1+(N_block-2)*floor((i-1)/N_block),(2*N_block-1))+1, :);
end

%compute Rc with the formula given in the master thesis
C_2 = kron(C, C);
Rc = C_2*T;

%compute rxhat using LS
rxhat = (Rc'*Rc)\Rc'*VecRy;

%the autocorrelation matrix Rx can be computed by the rxhat
VecRx = T*rxhat;
Rx = zeros(N_block);
for i = 1:N_block
   Rx(:,i) = VecRx((i-1)*N_block+1:i*N_block); 
end

%1. compute the spectrum from the dft of the autocorrelation rxhat(k)
Psp = fft(rxhat);
omega = 0:2*pi/length(rxhat):2*pi;
omega = omega(1:length(omega)-1);
figure()
plot(omega, abs(Psp))
hold on

%2. use minimum variance spectral estimate
e = NaN(length(Rx), 1);
Pmv = NaN(length(omega), 1);
for c1=1:length(omega)
    e(1) = 1;
for c2=2:length(Rx)
    e(c2) = exp(1i*(c2-1)*omega(c1));
end
denom = e'/Rx*e;
denom = real(denom);
Pmv(c1) = length(Rx)/denom;
end
plot(omega, Pmv)
set(gca,'XTick',[0 pi/2 pi 3*pi/2 2*pi])
set(gca,'XTickLabel',char('0', '\pi/2', '\pi', '3\pi/2','2\pi'))
xlim([0 2*pi])
xlabel('Radial frequency \omega')
ylabel('P_y(e^j^\omega)')
title('Power Spectrum Estimation Based on the Compressed Data')

%3. AR spectrum estimation
%create matrices for Yule-Walker equations
RxYW = Rx(2:end, 2:end);
rxYW = Rx(2:end, 1);

%run simulations for AR(1) up to AR(50)
maxP = 83;
AIC = zeros(maxP, 1);
%figure()
for p = 1:maxP
    RxP = RxYW(1:p, 1:p);
    rxP = rxYW(1:p);
    ap = -RxP\rxP;
    ep = RxP(1,1);
%     sum = 0;
    for k=1:p
       ep = ep + ap(k) * conj(rxP(k)); 
%        sum = sum + ap(k) * exp(-1i*k*omega);
    end
    AIC(p) = length(y_sig) * log10(ep) + 2 * p;
%     b0 = sqrt(ep);
%     if (mod(p,6) == 0)
%         %change figure after 6 plots, so that they do not get overcrowded
%         figure()
%     end
%     Par = ep./abs(1+sum);
%     plot(omega, abs(Par))
%     hold on
end

%plot the power spectrum estimated for an arbitrary order p = 10
p = 10;
RxP = RxYW(1:p, 1:p);
rxP = rxYW(1:p);
ap = -RxP\rxP;
ep = RxP(1,1);
sum = 0;
for k=1:p
   ep = ep + ap(k) * conj(rxP(k)); 
   sum = sum + ap(k) * exp(-1i*k*omega);
end
b0 = sqrt(ep);
Par = ep./abs(1+sum);
%figure()
plot(omega, abs(Par))
%hold on

%find the order p that minimizes the AIC and plot for this AR order
[M, I] = min(AIC);
p = I;
RxP = RxYW(1:p, 1:p);
rxP = rxYW(1:p);
ap = -RxP\rxP;
ep = RxP(1,1);
sum = 0;
for k=1:p
   ep = ep + ap(k) * conj(rxP(k)); 
   sum = sum + ap(k) * exp(-1i*k*omega);
end
b0 = sqrt(ep);
Par = ep./abs(1+sum);
plot(omega, abs(Par))
legend('Direct FFT', 'Minimum Variance', 'AR(10)', 'AR(56)')