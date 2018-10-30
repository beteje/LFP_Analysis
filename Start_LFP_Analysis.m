% Start File for LFP_Analysis
% set which analyses to perform
% then set the corresponding start parameters
% before calling LFP_Analysis

%% Which Analysis to Perform (1=yes 0=no)
PhaseLocking            =   1;
SpikeFieldCoherence     =   1;
OffsetPhaseLocking      =   0;
AmplitudeXCorr          =   0;
%% Spike (or LFP) to LFP (=1) or Spike (or LFP) to Multiple LFP (=0)
LFPSpike                =   0;

%% START PARAMETERS
DataFile        =   'FirstData.mat';    % Name of file data is stored in
FileName        =   'TestStSpike.mat';  % Name of file to store the results in
NumChannels     =   32;                 % Number of Channels
LFPRate         =   1e3;                % [Hz]
FiltFreq        =   300;                % Frequency for LPF [Hz]
MinNumSpikes    =   50;                 % Minimum number of spikes required for analysis
StandNumSpikes  =   1;                  % Set to 1 to Standardize the number of spikes to the minimum number of spikes
AnalyRange      =   [0;                 % Start times for block analysis [s]
                     120];              % End times for block analysis [s]
% Parameters for phase locking
PL_WavFreq      =   2.^((6:2:60)/8);    % Frequencies at which to perform phase locking analysis [Hz]
% Parameters for coherence
SFC_SegLength   =   480;                % Length of segments either side of spike for coherence analysis [ms]
% Parameters for offset phase locking (only used if LFPSpike=0)
OPL_FreqBand    =   [4,10];             % Minimum & Maximum frequencies for the frequency band used for the offset phase locking [Hz]
OPL_MaxLag      =   700;                % Maximum offset time to be analysed [ms]
OPL_StepSize    =   10;                 % Time steps to analyse the offsets at [ms]
% Parameters for Amplitude Cross Correlation
XCorr_MaxLag    =   100;                % Maximum lag for cross correlation [ms]
XCorr_FreqBands =   [7, 30;             % Start frequencies of frequency bands [Hz]
                     12,80];            % End frequencies of frequency bands [Hz]
XCorr_AnalyRange=   [0;                 % Start times of blocks for analysis [s]
                     120];              % End times of blocks for analysis [s]
XCorr_Channels  =   [1;                 % Pairs of channels to calculate the cross correlation between 
                     17];               % (only used if LFPSpike=1) otherwise all channels used

%% Call Analysis Program
Parameters  =   struct('Data_File',DataFile,...
                    'No_of_Channels',NumChannels,...
                    'LFP_Sampling_Frequency',LFPRate,...
                    'Low_Pass_Filter_Frequency',FiltFreq,...
                    'Frequencies_for_Phase_Locking_Analysis',PL_WavFreq,...
                    'Minimum_Required_Spikes',MinNumSpikes,...
                    'Standard_No_Spikes',StandNumSpikes',...
                    'Length_of_Segments_for_Coherence',SFC_SegLength,...
                    'AnalysisRange',AnalyRange,...
                    'Offset_Phase_Locking_Frequency_Band',OPL_FreqBand,...
                    'Maximum_Phase_Locking_Offset',OPL_MaxLag,...
                    'Phase_Locking_Offset_Step_Size',OPL_StepSize,...
                    'Amplitude_XCorr_Maximum_Lag',XCorr_MaxLag,...
                    'Amplitude_XCorr_Frequency_Bands',XCorr_FreqBands,...
                    'Amplitude_XCorr_Analysis_Range',XCorr_AnalyRange,...
                    'Amplitude_XCorr_Channel_Pairs',XCorr_Channels);

if LFPSpike
    LFP_Analysis(PhaseLocking,SpikeFieldCoherence,AmplitudeXCorr,FileName,Parameters);
else
    LFP_Analysis_XChannel(PhaseLocking,SpikeFieldCoherence,OffsetPhaseLocking,AmplitudeXCorr,FileName,Parameters);
end
