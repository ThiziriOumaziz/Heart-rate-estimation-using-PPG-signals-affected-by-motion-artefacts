function [f_sig] = BP(sig, Fs)

% Band-pass filter 
% sig : signal to filter
% Fs : sampling rate
% f_sig : Band-passed signal

    [b,a] = butter(2, [0.4 4]/(Fs/2),'bandpass');
    f_sig = filter(b,a,sig);
    
end
