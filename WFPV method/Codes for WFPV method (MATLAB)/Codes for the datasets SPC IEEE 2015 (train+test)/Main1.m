tic;
flops(0);
% Dataset IDs
IDData = {'DATA_01_TYPE01','DATA_02_TYPE02','DATA_03_TYPE02','DATA_04_TYPE02',...
    'DATA_05_TYPE02','DATA_06_TYPE02','DATA_07_TYPE02','DATA_08_TYPE02','DATA_09_TYPE02',...
    'DATA_10_TYPE02','DATA_11_TYPE02','DATA_12_TYPE02',...
    'TEST_S01_T01', 'TEST_S02_T01', 'TEST_S02_T02', 'TEST_S03_T02', ...
    'TEST_S04_T02', 'TEST_S05_T02', 'TEST_S06_T01', 'TEST_S06_T02',...
    'TEST_S07_T02', 'TEST_S08_T01'};

% Overall parameters
srate = 125;    % original sampling rate
FFTres = 1024;   % FFT resolution - user tunable parm1
WFlength = 15; % Wiener filter length - - user tunable parm2
allD = size(IDData,2); % number of recordings

% Filter parameters
CutoffFreqHzHP = 1; % 60 BPM;
CutoffFreqHzLP = 3;% 180 BPM
[b,a] = butter(4, [0.4 4]/(125/2),'bandpass');

% metrics
fullBPM=[];fullBPM0=[]; % for computing the correlation
myError = zeros(1,allD); % Beat per Minute errors
myErrorStd = zeros(1,allD); % std error
myRelError = zeros(1,allD); % relative error

% framework rule, 8s window 2s shift, no look into future
window   = 8 * srate;  % window length is 8 seconds
step     = 2 * srate;  % step size is 2 seconds


% loop for each recording
for idnb = 13 : 22
       
    % load the data
    load([IDData{idnb}]);
    if idnb>12
        ch1 = 1; ch2 = 2; ch3 = 3; ch4 = 4; ch5 = 5;
    else
        ch1 = 2; ch2 = 3; ch3 = 4; ch4 = 5; ch5 = 6;
    end
    windowNb = floor((length(sig)-window)/step) + 1;  % total number of windows(estimates)
    
    % initialization of variables
    BPM_est = zeros(1,windowNb); % estimated BPM
    rangeIdx = []; % range of search for the next estimates
    clear W1_FFTi W11_FFTi W2_FFTi W21_FFTi W1_PPG_ave_FFT_Clean W2_PPG_ave_FFT_Clean W11_PPG_ave_FFT_Clean PPG_ave_FFT_FIN W21_PPG_ave_FFT_Clean PPG_ave_FFT_FIN;
    
    
    for i =  [1 :  windowNb]
        curSegment = (i-1)*step+1 : (i-1)*step+window;
        curData = sig(:,curSegment);
        
        PPG1 = curData(ch1,:); PPG2 = curData(ch2,:);
        ACC_X = curData(ch3,:); ACC_Y = curData(ch4,:); ACC_Z = curData(ch5,:);
        
        % filtering
        PPG1 = filter(b,a,PPG1);
        PPG2 = filter(b,a,PPG2);
        ACC_X = filter(b,a,ACC_X);
        ACC_Y = filter(b,a,ACC_Y);
        ACC_Z = filter(b,a,ACC_Z);
        PPG_ave = 0.5*(PPG1-mean(PPG1))/(std(PPG1))+0.5*(PPG2-mean(PPG2))/(std(PPG2)); % mean overall
        
        % downsampling to 25Hz
        PPG_ave = downsample(PPG_ave,5);
        ACC_X = downsample(ACC_X,5);
        ACC_Y = downsample(ACC_Y,5);
        ACC_Z = downsample(ACC_Z,5);
        srate = 25; % new sampling rate
        
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
        
        Noise_FFT_norm = 1/3*(W1_ACC_X_FFT_norm+W1_ACC_Y_FFT_norm+W1_ACC_Z_FFT_norm);
        
               
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
                
        
        PPG_ave_FFT_FIN(i,:) = W1_PPG_ave_FFT_Clean(i,:)+ W2_PPG_ave_FFT_Clean(i,:) ; % ensambling W1 & W2
        
        hist_int = 25; % We start with +- 25 BPM history tracking window
        % History tracking
        if idnb>12
            % use the estimate of the max_abs_diff after 15 windows ~ 30s
            if i>15, hist_int = max(abs(diff(BPM_est(1:i-1))))+5; end; 
        else
            % use the estimate of the max_abs_diff after 30 windows ~ 60s
            if i>30, hist_int = max(abs(diff(BPM_est(1:i-1))))+5; end;        
        end
        
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
        if i>5 && abs(BPM_est(i)-BPM_est(i-1))>5
            %BPM_est(i) = BPM_est(i-1)+sign(BPM_est(i)-BPM_est(i-1))*10;
            ddd= polyfit(1:length(BPM_est(max(1,i-5):i-1)),BPM_est(max(1,i-5):i-1),1);
            BPM_est(i) = 0.8*BPM_est(i)+0.2*polyval(ddd,length(BPM_est(max(1,i-5):i-1))+1);
        end
        
        %Correction for delay
        mul=0.1;
        BPM_est(i) = BPM_est(i)+sum(sign(BPM_est(max(2,i-6):i)-BPM_est(max(1,i-7):i-1))*mul);
    end
    
    % Code to compute the error
    if idnb>12
        load(['True' IDData{idnb}(5:end)]);           % load groundtruth
    else
        load([IDData{idnb} '_BPMtrace']);           % load groundtruth
    end
    myError(idnb) = mean(abs(BPM0 - BPM_est(1:1:end)'));
    myRelError(idnb) = mean(abs(BPM0 - BPM_est(1:1:end)')./BPM0);
    myErrorStd(idnb) = std(abs(BPM0 - BPM_est(1:1:end)'));
    save(['Result_' IDData{idnb} '_DATA_WFPV'],'BPM_est', 'BPM0'); 
    
         
    if idnb > 0 fullBPM=[fullBPM BPM_est]; fullBPM0=[fullBPM0 BPM0']; end
    %plot best and worst performing patients
    %if idnb==13,          figure; plot(BPM0,'ro');hold on; plot(BPM_est(1:1:end),'o'); title(['Recording ' num2str(idnb)]); xlabel('Time (in frames)'); ylabel('HR'); legend({'Ground truth', 'Estimates'}); end;
    %if idnb==9,           figure; plot(BPM0,'ro');hold on; plot(BPM_est(1:1:end),'o'); title(['Recording ' num2str(idnb)]); xlabel('Time (in frames)'); ylabel('HR'); legend({'Ground truth', 'Estimates'}); end;
        
end
fprintf('Err12=%2.2f(%2.2f), Err11=%2.2f(%2.2f), ErrAll=%2.2f(%2.2f)', mean(myError(1:12)),mean(myErrorStd(1:12)),mean(myError(13:end)),mean(myErrorStd(13:end)),mean(myError),mean(myErrorStd));
for s=1:allD
    fprintf(' %2.2f', myError(s));
end

estimatedFLOPs = flops;
disp(['Estimated FLOPs: ' num2str(estimatedFLOPs)]);

elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);
% other plots
%BlandAltman(fullBPM0',fullBPM',{'Ground truth HR','Estimated HR'})
%tmp = corrcoef(fullBPM0,fullBPM); % correlation coefficient, offdiagonal
%hist(fullBPM0,20)
