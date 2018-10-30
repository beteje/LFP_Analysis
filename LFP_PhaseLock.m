function [SpikePhase,PhaseMean,MeanLength,PVal,PhaseLock] = LFP_PhaseLock(LFPData,SpikeLoc,WavFreq,Wav)
%% Calculate LFP phase for each spike using each of the supplied wavelets
%
%   Inputs
%   LFPData     --  local field potential data for a single channel
%   SpikeLoc    --  corresponding locations of the spikes
%   WavFreq     --  frequencies analysis is performed at
%   Wav         --  cell of wavelet coefficients
%
%   Outputs
%   SpikePhase  --  cell of phases for each spike at each frequency
%   PhaseMean   --  mean phase for each frequency
%   PVal        --  p values of each frequency
%   PhaseLock   --  frequencies at which the phase is locked to the spike

DataLength          =   length(LFPData);
NumFreq             =   length(Wav);
SpikePhase          =   cell(1,NumFreq); 
PhaseMean           =   zeros(1,NumFreq);
MeanLength          =   zeros(1,NumFreq);
PVal                =   zeros(1,NumFreq);
PhaseLock           =   [WavFreq;zeros(1,NumFreq)];
for j = 1:NumFreq
    x = conv(Wav{j},LFPData); % Compute wavelet transform
    x = x((1:DataLength) + (length(Wav{j})-1)/2); % Shift to remove delay
    
    AnalyticPhase = atan2(imag(x),real(x));
    
    % Find phase at each of the spike locations
    SpikePhase{j} = AnalyticPhase(logical(SpikeLoc));
    
    % Find phase mean for current frequency
    % compute weighted sum of cos and sin of angles
    r = sum(exp(1i*SpikePhase{j}));
    n = length(SpikePhase{j});
    % obtain mean angle by
    PhaseMean(j) = angle(r);
    % obtain mean length 
    MeanLength(j) = abs(r)./sum(ones(size(SpikePhase{j})));
    
    % Compute Rayleigh test for non-uniformity of circular data
    % compute Rayleigh's R 
    R = n*MeanLength(j);
    % compute p value using approximation 
    PVal(j) = exp(sqrt(1+4*n+4*(n^2-R^2))-(1+2*n));
    
    if PVal(j)<0.05/NumFreq
        PhaseLock(2,j) = 1;
    end
end