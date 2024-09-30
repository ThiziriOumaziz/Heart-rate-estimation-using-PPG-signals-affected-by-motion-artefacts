function curLoc = findClosestPeak(prevLoc, locs)
% Finds the closest peak to prevLoc

% prevLoc : frequency bin of previously located HR
% locs : frequency bins of peaks in search window

    % calculate distance to each peak
    dif = abs(prevLoc - locs);

    % find the closest peak
    [~, index] = min(dif);

    % return bin of closest peak
    curLoc = locs(index);

end