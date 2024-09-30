function curLoc = R1_findClosestPeak(prevLoc, locs, BPM_est, sauvv)
% Finds the closest peak to prevLoc

% prevLoc : frequency bin of previously located HR
% locs : frequency bins of peaks in search window

    % calculate distance to each peak
    dif = abs(prevLoc - locs);
    
    if (length(locs)>1)
    
        % find the closest peaks
        [difs, indexx] = mink(dif,2);
        locs = locs(indexx); 

        if difs(1)== difs(2) 

            %Observe the trend
            BPM_estimee = BPM_est(max(1,sauvv-10):sauvv-1);            
            BPMhatt = transpose(Smoother(BPM_estimee,5));

            dd= polyfit(1:length(BPMhatt),BPMhatt,1);
            deppBPM = polyval(dd,length(BPMhatt)+1);
            deppLoc = floor(((1024*(deppBPM))/(60*25))+1);

            [~,indd] = min(abs(locs-deppLoc));
            curLoc = locs(indd);

        else 
            curLoc = locs(1);
        end    
    else 
        % find the closest peak
        [~, index] = min(dif);

        % return bin of closest peak
        curLoc = locs(index);
    end        
end   
   
    
    

    
   