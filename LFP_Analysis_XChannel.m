function LFP_Analysis_XChannel(PhaseLocking,SpikeFieldCoherence,OffsetPhaseLocking,AmplitudeXCorr,FileName,Parameters)
%% Runs the requested analyses with the supplied parameters as LFP_Analysis 
%  except runs phase locking & SFC for all spikes with all LFP channels
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
if PhaseLocking == 1 || SpikeFieldCoherence==1 || OffsetPhaseLocking==1
    % Load the first channel 
    S = load(Parameters.Data_File,'FP01');
    % Find the number of samples in the LFP
    DataLength = size(S.FP01,1);
    % Create a new variable for the LFP data
    LFPData = zeros(DataLength,Parameters.No_of_Channels);
    % Store the filtered version of the LFP
    LFPData(:,1) = filtfilt(b,a,S.FP01(:,1));
    clear S;
    % Repeat for the remaining channels
    for i = 2:Parameters.No_of_Channels
        % Load LFP data for current channel
        CurrentLFP  =   sprintf('FP%02d',i);
        S           =   load(Parameters.Data_File,CurrentLFP);
        % Filter data
        LFPData(:,i)    =   filtfilt(b,a,S.(CurrentLFP)(:,1));
        clear S;
    end
        
    for i = 1:Parameters.No_of_Channels
        
        SpikeLetter = 0; % Create variable name for spike data
        Spike = strcat('SPK', sprintf('%02d',i), char(SpikeLetter+'a'));
        
        SpikeLoc    =   zeros(DataLength,1);
        CurrentChan =   sprintf('Ch%02d',i);
        
        % Check if Spike data exists
        while ismember(Spike,{SPKs.name})
            fprintf('Current Spike data: %s \n',Spike);
            
            S = load(Parameters.Data_File,Spike); % Load spike time series for current channel
            SpikeLoc(round(S.(Spike)*Parameters.LFP_Sampling_Frequency)) = 1; % Convert Spike Locations for LFP
            clear S
            CurrentChanSp = strcat(CurrentChan, char(SpikeLetter+'a'));
            
            for j=1:size(Parameters.AnalysisRange,2)
                % Calculate current data range
                Range = Parameters.AnalysisRange(1,j)*Parameters.LFP_Sampling_Frequency + 1: Parameters.AnalysisRange(2,j)*Parameters.LFP_Sampling_Frequency;
                % Check if there are at least the minimum required number of spikes otherwise move on to next section
                fprintf('\t Checking Number of spikes for %d-%d sec ... ',Parameters.AnalysisRange(1,j),Parameters.AnalysisRange(2,j));
                noSpikes = nnz(SpikeLoc(Range));
                if noSpikes < Parameters.Minimum_Required_Spikes
                    fprintf('Not enough spikes for analysis (only %d). Moving on to next section\n',noSpikes);
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
                
                % Create variable names for current analysis range
                SP  = strcat('SpikePhase_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                PM  = strcat('PhaseMean_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                ML  = strcat('MeanLength_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                PV  = strcat('PVal_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                PL  = strcat('PhaseLock_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                SF  = strcat('SFC_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                TA  = strcat('STA_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                FT  = strcat('fSTA_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                TP  = strcat('STP_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                OSP = strcat('OffsetSpikePhase_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                OPM = strcat('OffsetPhaseMean_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                OML = strcat('OffsetMeanLength_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                OZ  = strcat('OffsetZ_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                OPV = strcat('OffsetPVal_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                OPL = strcat('OffsetPhaseLock_', num2str(Parameters.AnalysisRange(1,j)), '_', num2str(Parameters.AnalysisRange(2,j)));
                
                for k=1:Parameters.No_of_Channels
                    %% Perform Phase Locking Analysis if selected
                    if PhaseLocking == 1
                        
                        % Use wavelets to find instantaneous phase for each spike at each frequency
                        [SpikePhase,PhaseMean,MeanLength,PVal,PhaseLock] = LFP_PhaseLock(LFPData_Range(:,k),SpikeLoc_Range,Parameters.Frequencies_for_Phase_Locking_Analysis,Wav);
                        
                        % Store results for current section of current channel;
                        Locking.(CurrentChanSp).(SP)(k,:)    =   SpikePhase;
                        Locking.(CurrentChanSp).(PM)(k,:)    =   PhaseMean;
                        Locking.(CurrentChanSp).(ML)(k,:)    =   MeanLength;
                        Locking.(CurrentChanSp).(PV)(k,:)    =   PVal;
                        if k==1
                            Locking.(CurrentChanSp).(PL)    =   PhaseLock;
                        else
                            Locking.(CurrentChanSp).(PL)(k+1,:) = PhaseLock(2,:);
                        end
                    end
                    
                    %% Perform Spike Field Coherence Analsyis if selected
                    if SpikeFieldCoherence == 1
                        % Calculate Spike Field Coherence
                        [SFC,STA,fSTA,STP,f] = LFP_Coherence(LFPData_Range(:,k),SpikeLoc_Range,(Parameters.Length_of_Segments_for_Coherence/1000)*Parameters.LFP_Sampling_Frequency,Parameters.LFP_Sampling_Frequency);
                        
                        % Store results for current section of current channel
                        Coherence.(CurrentChanSp).(SF)(k,:)  =   SFC;
                        Coherence.(CurrentChanSp).(TA)(k,:)  =   STA;
                        Coherence.(CurrentChanSp).(FT)(k,:)  =   fSTA;
                        Coherence.(CurrentChanSp).(TP)(k,:)  =   STP;
                        if i==1 && k==1
                            Coherence.f =   f;
                        end
                    end
                    
                    %% Perform Temporal Offset Phase Locking Analysis if selected
                    if OffsetPhaseLocking == 1
                        % Calculate the temporal phase spiking relationship 
                        [SpikePhase,PhaseMean,MeanLength,Z,PVal,PhaseLock] = LFP_OffsetPhaseLock(LFPData_Range(:,k),SpikeLoc_Range,Parameters.Offset_Phase_Locking_Frequency_Band(1),Parameters.Offset_Phase_Locking_Frequency_Band(2),Parameters.LFP_Sampling_Frequency,Parameters.Maximum_Phase_Locking_Offset,Parameters.Phase_Locking_Offset_Step_Size);
                        
                        % Store results for current section of current channel;
                        OffsetLocking.(CurrentChanSp).(OSP)(k,:)    =   SpikePhase;
                        OffsetLocking.(CurrentChanSp).(OPM)(k,:)    =   PhaseMean;
                        OffsetLocking.(CurrentChanSp).(OML)(k,:)    =   MeanLength;
                        OffsetLocking.(CurrentChanSp).(OZ)(k,:)     =   Z;
                        OffsetLocking.(CurrentChanSp).(OPV)(k,:)    =   PVal;
                        if k==1
                            OffsetLocking.(CurrentChanSp).(OPL)  =   PhaseLock;
                        else
                            OffsetLocking.(CurrentChanSp).(OPL)(k+1,:) = PhaseLock(2,:);
                        end
                    end
                end
                
                if PhaseLocking
                    save(FileName,'Locking','-append');
                end
                if SpikeFieldCoherence
                    save(FileName,'Coherence','-append');
                end
                if OffsetPhaseLocking
                    save(FileName,'OffsetLocking','-append');
                end
            end
            
            %% Prepare for next spike for current channel
            % Reset spike locations to zero
            SpikeLoc = zeros(DataLength,1);
            
            % Create next spike variable name
            SpikeLetter=SpikeLetter+1;
            Spike = strcat('SPK', sprintf('%02d',i), char(SpikeLetter+'a'));
        end % repeats until there is no more spike data for current channel
    end
end

%% Perform Amplitude Cross Correlation Analysis if selected
if AmplitudeXCorr ==1
    % Load the first channel
    S = load(Parameters.Data_File,'FP01');
    % Find the number of samples in the LFP
    DataLength = size(S.FP01,1);
    % Create a new variable for the LFP data
    LFPData = zeros(DataLength,Parameters.No_of_Channels);
    % Store the the LFP
    LFPData(:,1) = S.FP01(:,1);
    clear S;
    % Repeat for the remaining channels
    for i = 2:Parameters.No_of_Channels
        % Load LFP data for current channel
        CurrentLFP  =   sprintf('FP%02d',i);
        S           =   load(Parameters.Data_File,CurrentLFP);
        % Store data
        LFPData(:,i) = S.(CurrentLFP)(:,1);
        clear S;
    end
    
    for i=1:size(Parameters.Amplitude_XCorr_Frequency_Bands,2)
        fprintf('Current frequency band for Amplitude Cross Correlation: %d Hz to %d Hz\n',Parameters.Amplitude_XCorr_Frequency_Bands(1,i),Parameters.Amplitude_XCorr_Frequency_Bands(2,i));
        FreqBand = Parameters.Amplitude_XCorr_Frequency_Bands(:,i);
        CurrentFreq = sprintf('Freq_%d_%d',Parameters.Amplitude_XCorr_Frequency_Bands(1,i),Parameters.Amplitude_XCorr_Frequency_Bands(2,i));

        for j=1:size(Parameters.Amplitude_XCorr_Analysis_Range,2)
            fprintf('\t Performing Amplitude Cross Correlation Analysis for %d sec to %d sec\n',Parameters.Amplitude_XCorr_Analysis_Range(1,j),Parameters.Amplitude_XCorr_Analysis_Range(2,j));
            % Calculate current data range
            Range = Parameters.Amplitude_XCorr_Analysis_Range(1,j)*Parameters.LFP_Sampling_Frequency + 1: Parameters.Amplitude_XCorr_Analysis_Range(2,j)*Parameters.LFP_Sampling_Frequency;
            
            % Create variable names for current analysis range
            AXC  = strcat('Amp_XCorr_', num2str(Parameters.Amplitude_XCorr_Analysis_Range(1,j)), '_', num2str(Parameters.Amplitude_XCorr_Analysis_Range(2,j)));
            MAXCL  = strcat('Max_Amp_XCorr_Lag_', num2str(Parameters.Amplitude_XCorr_Analysis_Range(1,j)), '_', num2str(Parameters.Amplitude_XCorr_Analysis_Range(2,j)));
            AXCS  = strcat('Amp_XCorr_Significant_', num2str(Parameters.Amplitude_XCorr_Analysis_Range(1,j)), '_', num2str(Parameters.Amplitude_XCorr_Analysis_Range(2,j)));
            
            MaxXCorrLag = zeros(Parameters.No_of_Channels/2,Parameters.No_of_Channels/2);
            AmpXCorrSig = MaxXCorrLag;
%             for k=1:Parameters.No_of_Channels/2
            for k=1:2
                CurrentCh = sprintf('Ch%02d',k);
                fprintf('\t\t Amplitude Cross Correlations for %s\n',CurrentCh);
                
                % Calculate amplitude cross correlation between current channel and all other channels over current data range
                AmpXCorr = zeros(2*(Parameters.LFP_Sampling_Frequency/1000)*Parameters.Amplitude_XCorr_Maximum_Lag+1,Parameters.No_of_Channels/2);
%                 for l=Parameters.No_of_Channels/2+1:Parameters.No_of_Channels
                for l=17:19
                    [AmpXCorr(:,l-Parameters.No_of_Channels/2),Lags,MaxXCorrLag(k,l-Parameters.No_of_Channels/2),AmpXCorrSig(k,l-Parameters.No_of_Channels/2)] = LFP_AmpXCorr(LFPData(Range,k),LFPData(Range,l),Parameters.LFP_Sampling_Frequency,FreqBand,Parameters.Amplitude_XCorr_Maximum_Lag);
                end
                Amplitude_XCorr.(CurrentFreq).(CurrentCh).(AXC)     = AmpXCorr;
            end
            Amplitude_XCorr.(CurrentFreq).(MAXCL)   = MaxXCorrLag;
            Amplitude_XCorr.(CurrentFreq).(AXCS)    = AmpXCorrSig;
            Amplitude_XCorr.Lags = Lags;
            save(FileName,'Amplitude_XCorr','-append');
        end
    end
end
