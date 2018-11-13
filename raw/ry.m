function mean = ry( i, j, y_com, M_ruler )
    K = length(y_com)/M_ruler;
    sum = 0;
    for count=1:K
        sum = sum + y_com((count-1)*M_ruler+i+1)*conj(y_com((count-1)*M_ruler+j+1));
    end
    mean = sum/K;

end

