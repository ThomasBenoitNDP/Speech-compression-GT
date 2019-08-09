function [out_distance] = log_spectral_dist(sig1,sig2,L)
%Computing the log spectral distance between two signals

F1 =fft(sig1);   % F1(w)
F2 =fft(sig2);   % F2(w)

x1_cplx_cep = abs(ifft(log(F1)/log(10)));    
x2_cplx_cep = abs(ifft(log(F2)/log(10)));

c1 = x1_cplx_cep ;
c2 = x2_cplx_cep ;

distance = 2*sum((c1(1:L) - c2(1:L)).^2);
out_distance=distance;
end

