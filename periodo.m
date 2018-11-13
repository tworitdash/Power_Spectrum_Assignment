function perx = periodo(x, N)
  
        x_fft = fft(x);
        perx = (1/(2 * pi * N)) * abs(x_fft) .^ 2
    
end
