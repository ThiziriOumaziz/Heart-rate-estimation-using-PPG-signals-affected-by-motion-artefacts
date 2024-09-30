tic;
% Dataset IDs
IDData = {'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12','S13','S14','S15'};

% Overall parameters
srate_PPG = 64;  % original sampling rate of PPG
srate_ACC = 32; % original sampling rate of ACC
FFTres = 1024;   % FFT resolution - user tunable parm1
WFlength = 15; % Wiener filter length - - user tunable parm2
allD = size(IDData,2); % number of recordings

% Filter parameters
CutoffFreqHzHP = 0.65; % 39 BPM
CutoffFreqHzLP = 3.2; % 192 BPM
[b1,a1] = butter(4, [0.4 4]/(64/2),'bandpass'); %PPG
[b2,a2] = butter(4, [0.4 4]/(32/2),'bandpass'); %ACC

% metrics
fullBPM=[];fullBPM0=[]; % for computing the correlation
myError = zeros(1,allD); % Beat per Minute errors
myErrorStd = zeros(1,allD); % std error
myRelError = zeros(1,allD); % relative error

% framework rule, 8s window 2s shift, no look into future
window_PPG   = 8 * srate_PPG;  % window length is 8 seconds
window_ACC   = 8 * srate_ACC;
step_PPG = 2 * srate_PPG;  % step size is 2 seconds
step_ACC = 2 * srate_ACC;

% loop for each recording
for idnb = 15 : 15       
   % load the data
    load(['Dalia Data/' IDData{idnb}]);
    
    PPG_ = (pickle_data.signal.wrist.BVP (:,1)).';
    ACC_X_ = (pickle_data.signal.wrist.ACC (:,1)).';  
    ACC_Y_ = (pickle_data.signal.wrist.ACC (:,2)).';  
    ACC_Z_ = (pickle_data.signal.wrist.ACC (:,3)).'; 
    BPM0 = (pickle_data.label(1,:)).'; 
    
    windowNb = floor((length(PPG_)-window_PPG)/step_PPG) + 1;  % total number of windows(estimates)
    
        
    % initialization of variables
    BPM_est = zeros(1,windowNb); % estimated BPM
    rangeIdx = []; % range of search for the next estimates
    clear W1_FFTi W2_FFTi W1_PPG_ave_FFT_Clean W2_PPG_ave_FFT_Clean PPG_ave_FFT_FIN PPG_ave_FFT_FIN;
    
    
    for i =  [1 :  windowNb]
        curSegment_PPG = (i-1)*step_PPG+1 : (i-1)*step_PPG+window_PPG;
        curSegment_ACC = (i-1)*step_ACC+1 : (i-1)*step_ACC+window_ACC;
        PPG = PPG_(:,curSegment_PPG);
        ACC_X = ACC_X_(:,curSegment_ACC);
        ACC_Y = ACC_Y_(:,curSegment_ACC);
        ACC_Z = ACC_Z_(:,curSegment_ACC);
                
        % filtering
        PPG = filter(b1,a1,PPG);
        ACC_X = filter(b2,a2,ACC_X);
        ACC_Y = filter(b2,a2,ACC_Y);
        ACC_Z = filter(b2,a2,ACC_Z);
        PPG_ave =(PPG-mean(PPG))/(std(PPG)); % mean overall
        
        % downsampling to 32Hz
        PPG_ave = downsample(PPG_ave,2);
        srate = 32; % new sampling rate
        
        % Periodogram
        PPG_ave_FFT = fft(PPG_ave,FFTres);
        FreqRange = linspace(0,srate,size(PPG_ave_FFT,2));
        
        % finding the indices for the range of interest
        [~,lowR] = (min(abs(FreqRange-CutoffFreqHzHP)));
        [~,highR] = (min(abs(FreqRange-CutoffFreqHzLP)));
        
        %  Getting rid of most spectra outside the range of interest
        FreqRange = FreqRange(lowR:highR);
        PPG_ave_FFT = PPG_ave_FFT(lowR:highR);
        ACC_X_FFT= fft(ACC_X,FFTres); ACC_X_FFT = ACC_X_FFT(lowR:highR);
        ACC_Y_FFT= fft(ACC_Y,FFTres); ACC_Y_FFT = ACC_Y_FFT(lowR:highR);
        ACC_Z_FFT= fft(ACC_Z,FFTres); ACC_Z_FFT = ACC_Z_FFT(lowR:highR);
        
        
        % phase vocoder to refine spectral estimations
        FreqRangePPG = FreqRange;
        if i>1 % start phase vocoder for current and previous frames
            for ii=1:size(FreqRangePPG,2)
                curPhase = angle(PPG_ave_FFT(ii));
                prevPhase = angle(PPG_ave_FFTpr(ii)); vocoder = zeros(1,20);
                for n = 1:20
                    vocoder(n) = ((curPhase-prevPhase)+(2*pi*(n-1)))/(2*pi*2);
                end
                difference = vocoder - FreqRange(ii);
                [~, deltaidx] = min(abs(difference));
                FreqRangePPG(ii) = vocoder(deltaidx);
            end
        end
        
        % smooth phase vocoder frequency estimates
        FreqRangePPG = moving(FreqRangePPG,3);
        
        % save previous spectrum for the next phase vocoder call
        PPG_ave_FFTpr = PPG_ave_FFT;
        
        
        % Wiener filtering PPG-ACC, two types
             
        WC1 = WFlength; WC2 = WFlength;
        
        %Wiener 1 / abs & normalised
        W1_FFTi(i,:) = (abs(PPG_ave_FFT))/max(abs(PPG_ave_FFT));
        if i==1, W1_PPG_ave_FFT_ALL = W1_FFTi(i,:); else W1_PPG_ave_FFT_ALL = mean(W1_FFTi(max(1,i-WC1):i,:),1); end
        W1_PPG_ave_FFT_ALL_norm = (W1_PPG_ave_FFT_ALL)/max(W1_PPG_ave_FFT_ALL);
        W1_ACC_X_FFT_norm = (abs(ACC_X_FFT))/max(abs(ACC_X_FFT));
        W1_ACC_Y_FFT_norm = (abs(ACC_Y_FFT))/max(abs(ACC_Y_FFT));
        W1_ACC_Z_FFT_norm = (abs(ACC_Z_FFT))/max(abs(ACC_Z_FFT));
        WF1 = (1 - 1/3*(W1_ACC_X_FFT_norm+W1_ACC_Y_FFT_norm+W1_ACC_Z_FFT_norm)./(W1_PPG_ave_FFT_ALL_norm)); 
        WF1 (WF1<0) = -1; % limit negative -inf to -1
        W1_PPG_ave_FFT_Clean(i,:) = abs(PPG_ave_FFT).*WF1;
        
                
        %Wiener 2, abs & normalised
        W2_FFTi(i,:) = (abs(PPG_ave_FFT))/max(abs(PPG_ave_FFT));
        if i==1, W2_PPG_ave_FFT_ALL = W2_FFTi(i,:); else W2_PPG_ave_FFT_ALL = mean(W2_FFTi(max(1,i-WC2):i,:),1); end
        W2_PPG_ave_FFT_ALL_norm = (W2_PPG_ave_FFT_ALL)/max(W2_PPG_ave_FFT_ALL);
        W2_ACC_X_FFT_norm = (abs(ACC_X_FFT))/max(abs(ACC_X_FFT));
        W2_ACC_Y_FFT_norm = (abs(ACC_Y_FFT))/max(abs(ACC_Y_FFT));
        W2_ACC_Z_FFT_norm = (abs(ACC_Z_FFT))/max(abs(ACC_Z_FFT));
        WF2 = W2_PPG_ave_FFT_ALL_norm./(((W2_ACC_X_FFT_norm+W2_ACC_Y_FFT_norm+W2_ACC_Z_FFT_norm)/3)+W2_PPG_ave_FFT_ALL_norm); 
        W2_PPG_ave_FFT_Clean(i,:) = abs(PPG_ave_FFT).*WF2;
        W2_FFTi(i,:) = (W2_PPG_ave_FFT_Clean(i,:))/max(W2_PPG_ave_FFT_Clean(i,:));
        
               
        W1_PPG_ave_FFT_Clean(i,:) = W1_PPG_ave_FFT_Clean(i,:)/std(W1_PPG_ave_FFT_Clean(i,:)); 
        W2_PPG_ave_FFT_Clean(i,:) = W2_PPG_ave_FFT_Clean(i,:)/std(W2_PPG_ave_FFT_Clean(i,:)); 
       
    
        
        PPG_ave_FFT_FIN(i,:) = W1_PPG_ave_FFT_Clean(i,:)+ W2_PPG_ave_FFT_Clean(i,:); % ensambling W1 & W2
        
        hist_int = 25; % We start with +- 25 BPM history tracking window
        % History tracking
        
        % use the estimate of the max_abs_diff after 15 windows ~ 30s
        if i>15, hist_int = max(abs(diff(BPM_est(1:i-1))))+5; end 
        
        
        % HR estimation
        if isempty(rangeIdx)
            [extra, idx]= max(PPG_ave_FFT_FIN(i,:));
            BPM_est(i) = FreqRangePPG(idx(1))*60; 
            rangeIdx = idx(1)-round(hist_int/((FreqRange(2)-FreqRange(1))*60)):idx(1)+round(hist_int/((FreqRange(2)-FreqRange(1))*60));
        else
            [extra, idx]= max(PPG_ave_FFT_FIN(i,rangeIdx));
            BPM_est(i) = FreqRangePPG(rangeIdx(idx(1)))*60; 
            rangeIdx = rangeIdx(idx(1))-round(hist_int/((FreqRange(2)-FreqRange(1))*60)):rangeIdx(idx(1))+round(hist_int/((FreqRange(2)-FreqRange(1))*60));
        end
        rangeIdx(rangeIdx<1) = []; rangeIdx(rangeIdx>length(FreqRange)) = [];
        
        
        % Mild smoothing with linear regression prediction
        if i>5 && abs(BPM_est(i)-BPM_est(i-1))>10
            %BPM_est(i) = BPM_est(i-1)+sign(BPM_est(i)-BPM_est(i-1))*10;
            ddd= polyfit(1:length(BPM_est(max(1,i-5):i-1)),BPM_est(max(1,i-5):i-1),1);
            BPM_est(i) = 0.8*BPM_est(i)+0.2*polyval(ddd,length(BPM_est(max(1,i-5):i-1))+1);
        end
        
        %Correction for delay
        mul=0.1;
        BPM_est(i) = BPM_est(i)+sum(sign(BPM_est(max(2,i-6):i)-BPM_est(max(1,i-7):i-1))*mul);
    end
    
    % Code to compute the error
    
    myError(idnb) = mean(abs(BPM0 - BPM_est(1:1:end)'));
    myRelError(idnb) = mean(abs(BPM0 - BPM_est(1:1:end)')./BPM0);
    myErrorStd(idnb) = std(abs(BPM0 - BPM_est(1:1:end)'));
    save(['Result_' IDData{idnb} '_DALIA_WFPV'],'BPM_est', 'BPM0'); 
    
         
    if idnb > 0 fullBPM=[fullBPM BPM_est]; fullBPM0=[fullBPM0 BPM0']; end
    
    % plot best and worst performing patients    
    if idnb==5,  figure; plot(BPM0,'r');hold on; plot(BPM_est(1:1:end)); xlabel('Temps (en fenêtres)'); ylabel('Fréquence cardiaque (BPM)'); legend({'Vraie valeur', 'Valeur estimée'}); end;
    
end
fprintf('ErrAll=%2.2f(%2.2f)', mean(myError),mean(myErrorStd));
for s=1:allD
    fprintf(' %2.2f', myError(s));
end

toc 

% other plots
BlandAltman(fullBPM0',fullBPM',{'Ground truth HR','Estimated HR'})
tmp = corrcoef(fullBPM0,fullBPM); % correlation coefficient, offdiagonal
%hist(fullBPM0,20)
