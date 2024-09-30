function [pks,locations] = findPksInRange(clean_sig, R)
% Finds peaks of clean_sig in range R

% clean_sig : spectrum
% R : search range
% pks : peak values
% locations : peak frequency bins

    N = length(clean_sig);
    H = zeros(N,1);

    % filter for search range, handles cases when search range is negative
    H(R(R > 0)) = 1;
    % Band-pass filters spectrum to match search range
    clean_sig = clean_sig.*H;
    
    % returns peak values and locations
    [peaks,locations] = findpeaks(clean_sig); 
    [pks,ind] = maxk(peaks,3);
    locations = locations(ind); 
        
end