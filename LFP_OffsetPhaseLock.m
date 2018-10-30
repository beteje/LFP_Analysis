function [SpikePhase,PhaseMean,MeanLength,Z,PVal,PhaseLock] = LFP_OffsetPhaseLock(LFPData,SpikeLoc,f1,f2,Fs,MaxLag,StepSize)
%% Calculate LFP phase in one frequency band for the spikes with a range of temporal offsets
%
%   Inputs
%   LFPData     --  local field potential data for a single channel
%   SpikeLoc    --  corresponding locations of the spikes
%   f1,f2       --  minimum & maximum frequencies of the frequency band to be analysed
%   Fs          --  sampling frequency of the LFP data
%   MaxLag      --  maximum temporal offset
%   StepSize    --  step increments to analyse at
%
%   Outputs
%   SpikePhase  --  cell of phases for each spike at each offset
%   PhaseMean   --  mean phase for each offset
%   MeanLength  --  length of the mean phase for each offset
%   Z           --  Rayleighs Z score for each offset
%   PVal        --  p values of each offset
%   PhaseLock   --  offsets at which the phase is locked to the spike

StepSize    =   StepSize*(Fs/1000);
MaxLag      =   MaxLag*(Fs/1000);
NoLags      =   (MaxLag/StepSize)*2+1;
SpikePhase  =   cell(1,NoLags);
PhaseMean   =   zeros(1,NoLags);
MeanLength  =   zeros(1,NoLags);
Z           =   zeros(1,NoLags);
PVal        =   zeros(1,NoLags);
PhaseLock   =   [-MaxLag:StepSize:MaxLag;zeros(1,NoLags)];

b = fir1(2048,[f1,f2]/(Fs/2));
x = filtfilt(b,1,LFPData); 

AnalyticPhase = angle(hilbert(x));
SpikeLoc = [SpikeLoc,zeros(1,length(LFPData)-length(SpikeLoc))];

for i=1:NoLags
    % Shift the spikes by the desired offset
    if i<=MaxLag/StepSize
        tmpSpikeLoc = [zeros(MaxLag+(1-i)*StepSize,1);SpikeLoc(1:end-(MaxLag+(1-i)*StepSize))];
    elseif i>MaxLag/StepSize+1
        tmpSpikeLoc = [SpikeLoc((i-1)*StepSize-MaxLag+1:end);zeros((i-1)*StepSize-MaxLag,1)];
    else
        tmpSpikeLoc = SpikeLoc;
    end
    % Find phase at each of the adjusted spike locations
    SpikePhase{i} = AnalyticPhase(logical(tmpSpikeLoc));
    
    % Find phase mean for current frequency
    % compute weighted sum of cos and sin of angles
    r = sum(exp(1i*SpikePhase{i}));
    n = length(SpikePhase{i});
    % obtain mean angle by
    PhaseMean(i) = angle(r);
    % obtain mean length
    MeanLength(i) = abs(r)./sum(ones(size(SpikePhase{i})));
    
    % Compute Rayleigh test for non-uniformity of circular data
    % compute Rayleigh's R
    R = n*MeanLength(i);
    % compute Rayleigh's Z 
    Z(i) = R^2/n;

    % compute p value using approximation
    PVal(i) = exp(sqrt(1+4*n+4*(n^2-R^2))-(1+2*n));
    
    if PVal<0.05/NoLags
        PhaseLock(2,i) = 1;
    end
end