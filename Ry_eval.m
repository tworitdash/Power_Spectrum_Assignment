function ry_eval = Ry_eval(l, i, y_sig_r)
    load('lowpasssignal.mat');
    L = length(y_sig_r) / M_ruler;
    
    ry_eval_temp = 0; %temporary variable for correlation
    for m=1:L
        ry_eval_temp = ry_eval_temp + y_sig_r((m - 1) * M_ruler + l+ 1) * conj(y_sig_r((m - 1)*M_ruler + i +1));
    end
    %display(ry_eval_temp)
    ry_eval = ry_eval_temp / L;
end