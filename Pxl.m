function per_x = Pxl(x, win, n1, n2)
    load('lowpasssignal.mat')
    x = x(:);
    if nargin == 2
        n1 = 1; n2 = length(x); 
    end;
    per_x = 0;
    N = n2 - n1 + 1;
    L = floor(N / N_block);
    for i = 1:N_block-1
        w = ones(L, 1);
        if (win == 1) w = bartlett(L);
        elseif (win == 2) w = hamming(L);
        end;
        xw = x(n1:n1 + L-1).*w/norm(w);
        per_x = per_x + periodogram(xw);
        n1 = n1 + L;
    end
