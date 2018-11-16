function Px = bart(x,nsect,n_fft)
%
L = floor(length(x)/nsect);
Px = 0;
n1 = 1;
for i =1:nsect
    Px = Px + periodogram(x(n1:n1+L-1),n_fft)/nsect;
    n1 = n1 + L;
end
    