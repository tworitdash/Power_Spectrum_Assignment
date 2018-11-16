function perx = periodo(x, N)
  
        x_fft = fft(x);
        perx = (1/ N) * (x_fft) .* conj(x_fft);
    
end
