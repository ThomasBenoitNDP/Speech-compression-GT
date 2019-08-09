function [instantaneous_frequency] = freCh(analytic_signal,Fs)
    factor =  Fs/(2*pi);
    instantaneous_frequency = factor * diff(unwrap(angle(analytic_signal)));
    % Insert leading 0 in return-vector to maintain size
    instantaneous_frequency = [0 instantaneous_frequency];
end
