function per_x = Px(x, win, n1, n2)
%
    x = x(:);
    if nargin == 2
        n1 = 1; n2 = length(x); 
        display(n1, n2);
    end;
    N = n2 - n1 + 1
    w = ones(N, 1);
    if (win == 1) w = bartlett(N);
    elseif (win == 2) w = hamming(N);
    end;
    
    xw = x(n1:n2).*w/norm(w);
    per_x = N * periodogram(xw);
