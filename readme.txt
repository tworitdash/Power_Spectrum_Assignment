x_sig: Nyquist-rate version of the digital samples of the actual analog signal without the existence of white noise.

x_sig_and_noise:  Nyquist-rate version of the digital samples of the actual analog signal WITH the existence of additive white noise, 

stdev_noise: standard deviation of the noise.

The sub-Nyquist sampling in these files is simulated by implementing non-uniform selection (non-uniform sampling) on x_sig_and_noise. 


Applying this non-uniform sampling on one block of N_block consecutive samples leads to a block of M_ruler samples in y_sig. 
Hence we can say that the compression rate is given by M_ruler/N_block. 

ruler_index indicates which of N_block consecutive samples that are selected to form one block of M_ruler samples in y_sig.