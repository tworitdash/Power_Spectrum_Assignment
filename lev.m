function [a, epsilon] = lev(Rx_N_N, L)

    r = Rx_N_N(:);
    a = 1;
    p=L;
    epsilon=r(1);
    for j = 2:p+1
        gamma = -r(2:j)' * flipud(a)/epsilon;
        a = [a;0] + gamma*[0;conj(flipud(a))];
        epsilon=epsilon*(1 - abs(gamma)^2);
        
    end
end