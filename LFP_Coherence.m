function [SFC,STA,fSTA,STP,f] = LFP_Coherence(LFPData,SpikeLoc,NoPoints,SampFreq)
%% Calculate the Spike Field Coherence
%
%   Inputs
%   LFPData     --  local field potential data for a single channel
%   SpikeLoc    --  corresponding locations of the spikes
%   NoPoints    --  number of data points either side of the spike to be included in analysis
%   SampFreq    --  sampling frequency of LFP data
% 
%   Outputs
%   SFC         --  spike field coherence 
%   STA         --  spike triggered average
%   fSTA        --  frequency response of STA
%   STP         --  spike triggered power
%   f           --  frequencies the above was analysed at

SpikeInd = find(SpikeLoc); % Find indices of the Spike locations
SpikeInd = SpikeInd(SpikeInd>NoPoints & SpikeInd<length(LFPData)-NoPoints);

% Select segments of LFP data +- Segment Length around the Spike indices
Segments = LFPData(bsxfun(@plus, -NoPoints:NoPoints, SpikeInd));

% Calculate Spike Triggered Average
STA = mean(Segments);
% Calculate the frequency response of the STA
NFFT = 2^nextpow2(length(STA));
[fSTA,f] = pmtm(STA,4,NFFT,SampFreq);

% Calculate the frequency response of each spike segment
STP = zeros(size(fSTA));
for i = 1:length(SpikeInd);
    STP = STP + pmtm(Segments(i,:),4,NFFT,SampFreq);
end
% Calculate the Spike Triggered Power as a function of frequency
STP = STP/length(SpikeInd);
% Calculate the Spike Field Coherence as percentage
SFC = (fSTA./STP)*100;

