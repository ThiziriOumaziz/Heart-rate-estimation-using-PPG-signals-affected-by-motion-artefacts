% Dataset IDs
IDData = {'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12','S13','S14','S15'};

% Overall parameters
srate_PPG = 64;  % original sampling rate of PPG
srate_ACC = 32; % original sampling rate of ACC
N = 1024;   % FFT resolution - user tunable parm1
allD = size(IDData ,2); % number of recordings
lambda = 1e-10; %MMV 
MAX_ITERS = 4; % MMV : Number of iterations 
p = 0.8;
activity_type = 6; %activity type

% metrics
fullBPM=[];fullBPM0=[]; % for computing the correlation
myError = zeros(1,allD); % Beat per Minute errors
myErrorStd = zeros(1,allD); % std error
         
% framework rule, 8s window 2s shift, no look into future
window_PPG   = 8 * srate_PPG;  % window length is 8 seconds
window_ACC   = 8 * srate_ACC;
step_PPG = 2 * srate_PPG;  % step size is 2 seconds
step_ACC = 2 * srate_ACC;

% loop for each recording
for idnb = 6 : 6

    % load the data
    load(['Dalia Data/' IDData{idnb}]);
    
    activity  = pickle_data.activity(:,1);
    
          
    PPG_all = (pickle_data.signal.wrist.BVP (:,1)).';
    ACC_X_all = (pickle_data.signal.wrist.ACC (:,1)).';  
    ACC_Y_all = (pickle_data.signal.wrist.ACC (:,2)).';  
    ACC_Z_all = (pickle_data.signal.wrist.ACC (:,3)).'; 
    BPM0 = (pickle_data.label(1,:)).'; 
            
          
    %PPG indices 
    activity_64 = upsample (activity, 16);
    inter1 = find(activity_64==activity_type);
    PPG_init = min(inter1);
    PPG_end = max(inter1)+ 15; 
    
    %windows of interest
    
    w_Nb_init = floor(((PPG_init-1)-window_PPG)/step_PPG) + 2;
    w_Nb_end = floor(((PPG_end)-window_PPG)/step_PPG) + 1;
    
    
    windowNb = w_Nb_end - w_Nb_init + 1;
   
     
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
    
       
    init= w_Nb_init;
    for i =  w_Nb_init : w_Nb_end
        curSegment_PPG = (i-1)*step_PPG+1 : (i-1)*step_PPG+window_PPG;
        curSegment_ACC = (i-1)*step_ACC+1 : (i-1)*step_ACC+window_ACC;
        PPG_d = PPG_all(:,curSegment_PPG);
        ACC_X_d = ACC_X_all(:,curSegment_ACC);
        ACC_Y_d = ACC_Y_all(:,curSegment_ACC);
        ACC_Z_d = ACC_Z_all(:,curSegment_ACC);
                        
        % Filtering the series
        PPG_d = BP(PPG_d,srate_PPG);
        ACC_X = BP(ACC_X_d,srate_ACC);
        ACC_Y = BP(ACC_Y_d,srate_ACC);
        ACC_Z = BP(ACC_Z_d,srate_ACC);
              
        
        % downsampling to 32Hz
        PPG = downsample(PPG_d,2);
        srate=32;                  
                        
         
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
         [S, count] = MFOCUSS(Phi,Y,lambda,'p',p,'prune_gamma',1e-4,'max_iters',MAX_ITERS);
         S= abs(S).^2;
         
         %Quality of signal [0.8,2.5] Hz

        % finding the indices for the range of interest
        FreqRange = linspace(0,srate,N);
        [~,lowR] = (min(abs(FreqRange-0.8)));
        [~,highR] = (min(abs(FreqRange-2.5)));
        
        %Verification of the quality of the signal in the initialisation window 
        S_PPG = S(:,1);
        S_PPG = S_PPG(lowR:highR);
        Quality = kurtosis (S_PPG);
        Kurtosis (i-w_Nb_init+1) = Quality;
        
         %Spectral substraction
         [clean_sig] = Spect_sub(S, srate, N);
         
          
        if ((Quality < 10)&&(o==0))
            init=init+1;
            %mm = 1;
            continue;
        end   
         
         %SPT
                
         if i== init
             prevLoc = -1;
             prevBPM = -1;
         else 
             prevLoc = curLoc;
             prevBPM = curBPM;
         end  
         
         sauvv =i-w_Nb_init+1;                      
         [curLoc, curBPM, Trap_count,o,depLoc] = SPT_JOSS(clean_sig, srate, prevLoc, prevBPM, Trap_count, BPM_est,sauvv,init-w_Nb_init+1);
         BPM_est(i-w_Nb_init+1) = curBPM; 
    end
    
             
    % Code to compute the error
      
    myError(idnb) = mean(abs(BPM0(init:1:w_Nb_end) - BPM_est((init-w_Nb_init+1):1:end)'));
    myErrorStd(idnb) = std(abs(BPM0(init:1:w_Nb_end) - BPM_est((init-w_Nb_init+1):1:end)'));
        
    save(['Result_' IDData{idnb} '_DALIA_JOSS'],'BPM_est', 'BPM0');
    
    
   
end

fprintf('ErrAll=%2.2f(%2.2f)', mean(myError),mean(myErrorStd));
for s=1:allD
    fprintf(' %2.2f', myError(s));
end


    %if idnb > 0 fullBPM=[fullBPM BPM_est]; fullBPM0=[fullBPM0 BPM0']; end
    % plot best and worst performing patients
    %if idnb==14,           figure; plot(BPM0,'ro');hold on; plot(BPM_est(1:1:end),'o'); title(['Recording ' num2str(idnb)]); xlabel('Time (in frames)'); ylabel('HR'); legend({'Ground truth', 'Estimates'}); end;
    %if idnb==9,           figure; plot(BPM0,'ro');hold on; plot(BPM_est(1:1:end),'o'); title(['Recording ' num2str(idnb)]); xlabel('Time (in frames)'); ylabel('HR'); legen

         
         
       
        
    