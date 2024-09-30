function [curLoc, curBPM, Trap_count, o, depLoc] = SPT_JOSS(clean_sig, Fs, prevLoc, prevBPM, Trap_count, BPM_est, sauvv, init)
% Searches for the HR peak by using previously estimated HR values.

% clean_sig : cleaned spectrum
% Fs : sampling rate
% prevLoc : frequency bin of previously located HR, % initialized as -1
% prevBPM : previous HR, % initialized as -1
% Trap_count : counter that keeps track of the number of times % no peaks were found, initialized as 0
% curLoc : frequency bin of currently located HR
% curBPM : current HR

    % parameters SPT
    o = 1; 
    delta1 = 15;
    delta2 = 25;
    delta3 = 30;
    N = length(clean_sig);
    FreqRange = linspace(0,Fs,N);
    sauveur = zeros(1024,1);
    depLoc=0;

    % finding the indices for the range [0.8,2.5]
    [~,lowR] = (min(abs(FreqRange-0.8)));
    [~,highR] = (min(abs(FreqRange-2.5)));

     
    % initialization state
    if (prevLoc == -1 && prevBPM == -1) 
        sauveur(lowR:highR) = clean_sig(lowR:highR);
        [~,curLoc] = max(sauveur);
        curBPM = 60 * (curLoc - 1) / N * Fs;
    else
        
    % set search range
        R1 = (prevLoc - delta1) : (prevLoc + delta1);
        [~,locs] = findPksInRange(clean_sig, R1);
        numOfPks = length(locs);
        if (numOfPks >= 1)
        % find closest peak
            curLoc = R1_findClosestPeak(prevLoc, locs, BPM_est, sauvv);
            curBPM = 60 * (curLoc - 1) / N * Fs;
            if (abs(prevBPM - curBPM) > 12)
                curLoc = prevLoc;
                curBPM = prevBPM;
            end    
                
        else
            % increase search range
            R2 = (prevLoc - delta2) : (prevLoc + delta2);
            % find peaks in range2
            [pks,locs] = findPksInRange(clean_sig, R2);
            numOfPks = length(locs);
            if (numOfPks >= 1)
            % find maximum peak
                [~, maxIndex] = max(pks);
                curLoc = locs(maxIndex);
                curBPM = 60 * (curLoc - 1) / N * Fs;
                
                if (abs(prevBPM - curBPM) > 12)
                    curLoc = prevLoc;
                    curBPM = prevBPM;
                end   
                            
            else
            % choose prev BPM
                curLoc = prevLoc;
                curBPM = prevBPM;
            end
        end
    end
    
    % validation   
        
    if (curLoc == prevLoc)
        Trap_count = Trap_count + 1;
        if (Trap_count > 2)
            
        % discover
                                  
            BPM_estime = BPM_est(max(init,sauvv-30):sauvv-1);            
            BPMhat = transpose(Smoother(BPM_estime,20));
            
            ddd= polyfit(1:length(BPMhat),BPMhat,1);
            depBPM = polyval(ddd,length(BPMhat)+1);
             
            depLoc = floor(((1024*(depBPM))/(60*25))+1);

            R3 = depLoc-delta3:depLoc+delta3;
            
            [~,locs] = findPksInRange(clean_sig, R3);
            curLoc = findClosestPeak(prevLoc, locs);
            
           
            if (isempty(curLoc))
                curLoc = depLoc;
            end    
                      
            curBPM = 60 * (curLoc - 1) / N * Fs;
        end
    else
        Trap_count = 0;
    end
end

%depBPM = BPMhat(length(BPMhat));

%[pkks,loccs] = findPksInRange(clean_sig, R3);
                      
% find maximum peak
%[~, maxIndexx] = max(pkks);
%curLoc = loccs(maxIndexx);
                       
