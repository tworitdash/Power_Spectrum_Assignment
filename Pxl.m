function per_x = Pxl(x, N_block)

    X = zeros(length(N_block) ,1);
    for i = 1:fix(length(x)/N_block)-1
        
       x_i = x(i*N_block+1 : (i+1)*N_block);
       X_i = fft(x_i);
       X_i = abs(X_i) .^ 2;
       X = X + X_i;
    end
    per_x = X/length(x);
end