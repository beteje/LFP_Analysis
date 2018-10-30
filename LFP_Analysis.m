function LFP_Analysis(PhaseLocking,SpikeFieldCoherence,AmplitudeXCorr,FileName,Parameters)
%% Runs the requested analyses with the supplied parameters 
%   Inputs
%   PhaseLocking            --  whether or not to perform phase locking analysis
%   SpikeFieldCoherence     --  whether or not to perform spike field coherence analysis
%   FileName                --  save file for the results
%   Parameters              --  start parameters set in Start_LFP_Analysis
%
%   Outputs
%   all data is saved in FileName
%       Parameters                              --  initial start parameters
%   for channel number xx & time block from time1 to time 2
%       Locking.Chanxx.SpikePhase_time1_time2   --  cell of phases for each spike at each frequency 
%       Locking.Chanxx.PhaseMean_time1_time2    --  mean phase for each frequency
%       Locking.Chanxx.PhaseLength_time1_time2  --  mean length for each frequency
%       Locking.Chanxx.PVal_time1_time2         --  p values of each frequency 
%       Locking.Chanxx.PhaseLock_time1_time2    --  frequencies at which the phase is locked to the spike and corresponding phase for 
%       Coherence.Chanxx.SFC_time1_time2        --  spike field coherence 
%       Coherence.Chanxx.STA_time1_time2        --  spike triggered average
%       Coherence.Chanxx.fSTA_time1_time2       --  frequency response of STA
%       Coherence.Chanxx.STP_time1_time2        --  spike triggered power
%       Coherence.f                             --  frequencies the coherence was analysed at

% Create low pass filter ready for analysis 
[b,a]   =   butter(4,Parameters.Low_Pass_Filter_Frequency/Parameters.LFP_Sampling_Frequency);
if PhaseLocking 
    % Create the Morlets ready for analysis
    Wav     =   LFP_GenMorlets(Parameters.Frequencies_for_Phase_Locking_Analysis,Parameters.LFP_Sampling_Frequency);
end
% Find all the Spike data in data file
SPKs    =   whos('-file',Parameters.Data_File,'-regexp', 'SPK\w*');
% Save start parameters
if exist(FileName,'file')
    save(FileName,'Parameters','-append');
else
    save(FileName,'Parameters');
end

%% Load Data if either phase locking or spike field coherence is selected
if PhaseLocking == 1 || SpikeFieldCoherence==1
    for i = 1:Parameters.No_of_Channels
        
        k = 0; % Create variable name for spike data
        Spike = strcat('SPK', sprintf('%02d',i), char(k+'a'));
        
        % Load LFP data for current channel
        CurrentLFP  =   sprintf('FP%02d',i);
        S           =   load(Parameters.Data_File,CurrentLFP);
        LFPData     =   S.(CurrentLFP)(:,1);
        clear S;
        
        % Filter data
        LFPData     =   filtfilt(b,a,LFPData);
        
        DataLength  =   length(LFPData);
        SpikeLoc    =   zeros(DataLength,1);
        CurrentChan =   sprintf('Ch%02d',i);
        
        % Check if Spike data exists
        while ismember(Spike,{SPKs.name})
            fprintf('Current Spike data: %s \n',Spike);
            
            S = load(Parameters.Data_File,Spike); % Load spike time series for current channel
            SpikeLoc(round(S.(Spike)*Parameters.LFP_Sampling_Frequency)) = 1; % Convert Spike Locations for LFP
            clear S
            CurrentChanSp = strcat(CurrentChan, char(k+'a'));
            
            for j=1:size(Parameters.AnalysisRange,2)
                % Calculate current data range
                Range = Parameters.AnalysisRange(1,j)*Parameters.LFP_Sampling_Frequency + 1: Parameters.AnalysisRange(2,j)*Parameters.LFP_Sampling_Frequency;
                % Check if there are at least the minimum required number of spikes otherwise move on to next section
                fprintf('\t Checking Number of spikes for %d-%d sec ... ',Parameters.AnalysisRange(1,j),Parameters.AnalysisRange(2,j));
                noSpikes = nnz(SpikeLoc(Range));
                if noSpikes < Parameters.Minimum_Required_Spikes
                    fprintf('Not enough spikes for analysis (only %d). Moving on to next section\n',nnz(SpikeLoc(Range)));
                    continue;
                else
                    fprintf('OK\n');
                    LFPData_Range   = LFPData(Range,:);
                    SpikeLoc_Range  = SpikeLoc(Range);
                end
                
                 % Select a standardized number of spikes
                if Parameters.Standard_No_Spikes == 1
                    RemoveSpikes    = randperm(noSpikes,noSpikes-Parameters.Minimum_Required_Spikes);
                    SpikeInd        = find(SpikeLoc_Range);
                    SpikeLoc_Range(SpikeInd(RemoveSpikes)) = 0;
                end
                
                %% Perform Phase Locking Analysis if selected
                if PhaseLocking == 1
                    fprintf('\t\t Performing Phase Locking Analysis\n');
                    
                    % Use wavelets to find instantaneous phase for each spike at each frequency
                    [SpikePhase,PhaseMean,MeanLength,PVal,PhaseLock] = LFP_PhaseLock(LFPData_Range,SpikeLoc_Range,Parameters.Frequencies_for_Phase_Locking_Analysis,Wav);
                    
                    % Store results for current section of current channel;
                    SP = strcat('SpikePhase_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                    PM = strcat('PhaseMean_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                    ML = strcat('MeanLength_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                    PV = strcat('PVal_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                    PL = strcat('PhaseLock_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                    Locking.(CurrentChanSp).(SP)    =   SpikePhase;
                    Locking.(CurrentChanSp).(PM)    =   PhaseMean;
                    Locking.(CurrentChanSp).(ML)    =   MeanLength;
                    Locking.(CurrentChanSp).(PV)    =   PVal;
                    Locking.(CurrentChanSp).(PL)    =   PhaseLock;
                    save(FileName,'Locking','-append');
                end
                
                %% Perform Spike Field Coherence Analsyis if selected
                if SpikeFieldCoherence == 1
                    fprintf('\t\t Performing Spike Field Coherence Analysis\n');
                    % Calculate Spike Field Coherence
                    [SFC,STA,fSTA,STP,f] = LFP_Coherence(LFPData_Range,SpikeLoc_Range,(Parameters.Length_of_Segments_for_Coherence/1000)*Parameters.LFP_Sampling_Frequency,Parameters.LFP_Sampling_Frequency);
                    
                    % Store results for current section of current channel
                    SF = strcat('SFC_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                    TA = strcat('STA_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                    FT = strcat('fSTA_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                    TP = strcat('STP_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                    Coherence.(CurrentChanSp).(SF)  =   SFC;
                    Coherence.(CurrentChanSp).(TA)  =   STA;
                    Coherence.(CurrentChanSp).(FT)  =   fSTA;
                    Coherence.(CurrentChanSp).(TP)  =   STP;
                    Coherence.f                     =   f;
                    save(FileName,'Coherence','-append');
                end
            end
            
            %% Prepare for next spike for current channel
            % Reset spike locations to zero
            SpikeLoc = zeros(DataLength,1);
            
            % Create next spike variable name
            k=k+1;
            Spike = strcat('SPK', sprintf('%02d',i), char(k+'a'));
        end % repeats until there is no more spike data for current channel
    end
end

%% Perform Amplitude Cross Correlation Analysis if selected
if AmplitudeXCorr ==1
    for i=1:size(Parameters.Amplitude_XCorr_Frequency_Bands,2)
        fprintf('Current frequency band for Amplitude Cross Correlation: %d Hz to %d Hz\n',Parameters.Amplitude_XCorr_Frequency_Bands(1,i),Parameters.Amplitude_XCorr_Frequency_Bands(2,i));
        FreqBand = Parameters.Amplitude_XCorr_Frequency_Bands(:,i);
        CurrentFreq = sprintf('Freq_%d_%d',Parameters.Amplitude_XCorr_Frequency_Bands(1,i),Parameters.Amplitude_XCorr_Frequency_Bands(2,i));
    
        for j=1:size(Parameters.Amplitude_XCorr_Channel_Pairs,2)
            Ch = Parameters.Amplitude_XCorr_Channel_Pairs(:,j);
            betweenChans = sprintf('Ch%02dCh%02d',Ch(1),Ch(2));
            
            % Load first channel for amplitude cross correlation analysis 
            CurrentLFP  =   sprintf('FP%02d',Ch(1));
            S           =   load(Parameters.Data_File,CurrentLFP);
            LFPData1    =   S.(CurrentLFP)(:,1); 
            clear S;
            fprintf('\t Amplitude Cross Correlation between %s',CurrentLFP);
    
            % Load second channel for information flow analysis
            CurrentLFP  =   sprintf('FP%02d',Ch(2));
            S           =   load(Parameters.Data_File,CurrentLFP);
            LFPData2    =   S.(CurrentLFP)(:,1); 
            clear S;
            fprintf(' and %s\n',CurrentLFP);
            
            AmpXCorr = zeros(2*(Parameters.LFP_Sampling_Frequency/1000)*Parameters.Amplitude_XCorr_Maximum_Lag+1,size(Parameters.AnalysisRange,2));
            MaxXCorrLag = zeros(1,size(Parameters.Amplitude_XCorr_Analysis_Range,2));
            AmpXCorrSig = MaxXCorrLag;
            for k=1:size(Parameters.Amplitude_XCorr_Analysis_Range,2)
                fprintf('\t\t Performing Amplitude Cross Correlation Analysis for %d sec to %d sec\n',Parameters.Amplitude_XCorr_Analysis_Range(1,k),Parameters.Amplitude_XCorr_Analysis_Range(2,k));
                % Calculate current data range
                Range = Parameters.Amplitude_XCorr_Analysis_Range(1,k)*Parameters.LFP_Sampling_Frequency + 1: Parameters.Amplitude_XCorr_Analysis_Range(2,k)*Parameters.LFP_Sampling_Frequency;
                % Calculate the amplitude cross correlation for the current pair of channels over current data range
                [AmpXCorr(:,k),Lags,MaxXCorrLag(:,k),AmpXCorrSig(:,k)] = LFP_AmpXCorr(LFPData1(Range),LFPData2(Range),Parameters.LFP_Sampling_Frequency,FreqBand,Parameters.Amplitude_XCorr_Maximum_Lag);
            end
            
            Amplitude_XCorr.(CurrentFreq).(betweenChans).Amp_XCorr = AmpXCorr;
            Amplitude_XCorr.(CurrentFreq).(betweenChans).Max_Amp_XCorr_Lag = MaxXCorrLag;
            Amplitude_XCorr.Lags = Lags;
            Amplitude_XCorr.(CurrentFreq).(betweenChans).Amp_XCorr_Significant = AmpXCorrSig;
            save(FileName,'Amplitude_XCorr','-append');
        end  
    end
end
