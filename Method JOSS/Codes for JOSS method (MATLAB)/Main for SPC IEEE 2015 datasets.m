tic;
% Dataset IDs
IDData = {'DATA_01_TYPE01','DATA_02_TYPE02','DATA_03_TYPE02','DATA_04_TYPE02',...
    'DATA_05_TYPE02','DATA_06_TYPE02','DATA_07_TYPE02','DATA_08_TYPE02','DATA_09_TYPE02',...
    'DATA_10_TYPE02','DATA_11_TYPE02','DATA_12_TYPE02',...
    'TEST_S01_T01', 'TEST_S02_T01', 'TEST_S02_T02', 'TEST_S03_T02', ...
    'TEST_S04_T02', 'TEST_S05_T02', 'TEST_S06_T01', 'TEST_S06_T02',...
    'TEST_S07_T02', 'TEST_S08_T01'};

% Overall parameters
Fs = 125;    % original sampling rate
N = 1024;   % FFT resolution - user tunable parm1
allD = size(IDData,2); % number of recordings
lambda = 1e-10; %MMV 
MAX_ITERS = 3; % MMV : Number of iterations 


% metrics
fullBPM=[];fullBPM0=[]; % for computing the correlation
myError = zeros(1,allD); % Beat per Minute errors
myErrorStd = zeros(1,allD); % std error

% framework rule, 8s window 2s shift, no look into future
window   = 8 * Fs;  % window length is 8 seconds
step     = 2 * Fs;  % step size is 2 seconds

% loop for each recording
for idnb = 10 : 10
    % load the data
    load([IDData{idnb}]);
    if idnb>12
        ch1 = 1; ch2 = 2; ch3 = 3; ch4 = 4; ch5 = 5;
    else
        ch1 = 2; ch2 = 3; ch3 = 4; ch4 = 5; ch5 = 6;
    end
    
    Fs=125;
            
    PPG_all_1 = sig(ch1,:); PPG_all_2 = sig(ch2,:); ACC_X_all = sig(ch3,:); ACC_Y_all = sig(ch4,:); ACC_Z_all = sig(ch5,:);
      
    windowNb = floor((length(PPG_all_1)-window)/step) + 1;
    
    % initialization of variables
    BPM_est = zeros(1,windowNb); % estimated BPM
    Kurtosis = zeros(1,windowNb);
    clean_sig = zeros(1024,1);
    count = 0;
    o = 0;
    Trap_count = 0;
    curBPM = 0;
    curLoc = 0;
    depLoc = 0;
    prevLoc = 0;
    prevBPM =0;
        
    init=1;
    for i = [1 : windowNb]
        curSegment = (i-1)*step+1 : (i-1)*step+window;
       
        PPG_d_1 = PPG_all_1(1,curSegment); PPG_d_2 = PPG_all_2(1,curSegment); ACC_X_d = ACC_X_all(1, curSegment); ACC_Y_d = ACC_Y_all(1, curSegment); ACC_Z_d = ACC_Z_all(1, curSegment);
                        
        % Filtering the series
        Fs=125;
        PPG_d_1 = BP(PPG_d_1,Fs);
        PPG_d_2 = BP(PPG_d_2,Fs);
        ACC_X_d = BP(ACC_X_d,Fs);
        ACC_Y_d = BP(ACC_Y_d,Fs);
        ACC_Z_d = BP(ACC_Z_d,Fs);
        PPG_d = 0.5*(PPG_d_1-mean(PPG_d_1))/(std(PPG_d_1))+0.5*(PPG_d_2-mean(PPG_d_2))/(std(PPG_d_2));
        
        
        % downsampling to 25Hz
        PPG = downsample(PPG_d,5);
        ACC_X = downsample(ACC_X_d,5);
        ACC_Y = downsample(ACC_Y_d,5);
        ACC_Z = downsample(ACC_Z_d,5);
        Fs=25;                  
                        
         
        %MMV Model 
        
        % construct fourier matrix
        M = length(PPG);
        Phi = zeros(M,N);
        complex_factor = 1i*2*pi/N;
        for m = 1:M
            for n = 1:N
                Phi(m,n) = exp(complex_factor*(m-1)*(n-1));
            end
        end
        
        %construct Y 
         Y = zeros (M,4);
         Y(:,1) = PPG;
         Y(:,2) = ACC_X;
         Y(:,3) = ACC_Y;
         Y(:,4) = ACC_Z;
         
                  
         %Compute the joint sparse spectra
         [S, count] = MFOCUSS(Phi,Y,lambda,'p',0.9,'prune_gamma',1e-4,'max_iters',MAX_ITERS);
         S= abs(S).^2;
         
         %Quality of signal [0.8,2.5] Hz

        % finding the indices for the range of interest
        FreqRange = linspace(0,Fs,N);
        [~,lowR] = (min(abs(FreqRange-0.8)));
        [~,highR] = (min(abs(FreqRange-2.5)));
        
        %Verification of the quality of the signal in the initialisation window 
        S_PPG = S(:,1);
        S_PPG = S_PPG(lowR:highR);
        Quality = kurtosis (S_PPG);
        Kurtosis (i) = Quality;
        
         %Spectral substraction
         [clean_sig, S] = Spect_sub(S, Fs, N);
         
          
        if ((Quality < 10)&&(o==0))
            init=init+1;
            %mm = 1;
            continue;
        end   
         
         %SPT
                
         if i==init
             prevLoc = -1;
             prevBPM = -1;
         else 
             prevLoc = curLoc;
             prevBPM = curBPM;
         end  
         
         sauvv = i;                      
         [curLoc, curBPM, Trap_count,o,depLoc] = SPT_JOSS(clean_sig, Fs, prevLoc, prevBPM, Trap_count, BPM_est,sauvv,init);
         BPM_est(i) = curBPM; 
    end
    
             
    % Code to compute the error
    if idnb>12
        load(['True' IDData{idnb}(5:end)]);           % load groundtruth
    else
        load([IDData{idnb} '_BPMtrace']);           % load groundtruth
    end
    myError(idnb) = mean(abs(BPM0(init:1:end)- BPM_est(init:1:end)'));
    myErrorStd(idnb) = std(abs(BPM0(init:1:end) - BPM_est(init:1:end)'));
    save(['Result_' IDData{idnb} '_DATA_JOSS'],'BPM_est', 'BPM0');
    
end
    fprintf('Err1=%2.2f(%2.2f), Err2=%2.2f(%2.2f), ErrAll=%2.2f(%2.2f)', mean(myError(1:12)),mean(myErrorStd(1:12)),mean(myError(13:end)),mean(myErrorStd(13:end)),mean(myError),mean(myErrorStd));
    for s=1:allD
        fprintf(' %2.2f', myError(s));
    end


    if idnb > 0 fullBPM=[fullBPM BPM_est]; fullBPM0=[fullBPM0 BPM0']; end
    % plot best and worst performing patients
    %if idnb==13,           figure; plot(BPM0,'rx');hold on; plot(BPM_est(1:1:end),'x'); xlabel('Temps (en fenêtres)'); ylabel('Fréquence cardiaque (BPM)'); legend({'Vraie valeur', 'Valeur estimée'}); end;
    if idnb==22,           figure; plot(BPM0,'r');hold on; plot(BPM_est(1:1:end)); xlabel('Temps (en fenêtres)'); ylabel('Fréquence cardiaque (BPM)'); legend({'Vraie valeur', 'Valeur estimée'}); end;

elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);         
                
       
        
    