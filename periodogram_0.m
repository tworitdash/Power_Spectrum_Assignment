function Px = periodogram_0(x,n_fft,n1,n2)
%
    x = x(:);
    if nargin == 2 
        n1 = 1; n2 = length(x); end
    Px = abs(fft(x(n1:n2),n_fft)).^2/(n2-n1);
    Px(1) = Px(2);
end
