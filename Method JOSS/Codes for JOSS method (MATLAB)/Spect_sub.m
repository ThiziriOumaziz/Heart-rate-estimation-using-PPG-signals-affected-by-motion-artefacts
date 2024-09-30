function [clean_sig, S] = Spect_sub(S, Fs, N)
% Removes MA from PPG by means of spectral subtraction

% S : Matrix of Sparse Spectra PPG + 3 ACC
% Fs : sampling rate
% N : resolution
% clean_sig : cleansed sparse spectrum of PPG

    FreqRange = linspace(0,Fs,N);
    aggression = 0.99;

    clean_sig = S(:,1);
    clean_sig = clean_sig/max(clean_sig);
    
    for i = 1:4
        SS = S(:,i);
    % normalize spectra
        SS = SS/max(SS);
        S(:,i) = SS;
    end

    for i = 1:length(FreqRange)
    % max of acceleration at f_i
        C(i) = max([S(i,2),S(i,3),S(i,4)]);
    % spectral subtraction
        clean_sig(i) = clean_sig(i) - aggression*C(i);
    end

    % set values smaller than 1/4-th of the max to 0
    p_max = max(clean_sig);
    clean_sig(clean_sig < p_max/4) = 0;

end
