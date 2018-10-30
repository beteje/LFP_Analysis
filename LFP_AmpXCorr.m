function [AmpXCorr,Lags,MaxXCorrLag,AmpXCorrSig] = LFP_AmpXCorr(LFPData1,LFPData2,Fs,FreqBand,MaxLag)
%% Calculates the cross correlation of the amplitude envelope between two LFP channels in specified frequency band
% 
% Inputs
% LFPData 1&2   -- LFP data from different regions to be compared
% Fs            -- Sampling frequency of the LFP data
% FreqBand      -- Minium and maximum frequency of the frequency band
% MaxLag        -- Maximum lag to use in cross correlation
%
% Outputs
% AmpXCorr      -- Cross correlation values of the amplitude envelopes 
% Lags          -- Lags at which cross correlation was calculates
% MaxXCorrLag   -- Lag at which the cross correlation is maximum (negative value indicates LFPData1 leads LFPData2 and vice versa)

% Filter LFP data into desired frequency band                            
b   = fir1(2048,FreqBand/(Fs/2));
x1  = filtfilt(b,1,LFPData1);                    
x2  = filtfilt(b,1,LFPData2);                  

% Calculate the instantaneous amplitude of the filtered data (removing the mean)
InstAmp1        = abs(hilbert(x1));              
InstAmp1        = InstAmp1-mean(InstAmp1);                        

InstAmp2        = abs(hilbert(x2));                          
InstAmp2        = InstAmp2-mean(InstAmp2);

% Calculate the cross correlations and identify lag at which the cross correlation is maximum
MaxLag          = (Fs/1000)*MaxLag;
[AmpXCorr,Lags] = xcorr(InstAmp1,InstAmp2,MaxLag,'coeff');    
Lags            = (Lags./Fs)*1000; % converts lags to miliseconds
MaxXCorrLag     = Lags(AmpXCorr==max(AmpXCorr));     

% Calculate significance of result using bootstrap procedure
% Create 1000 random shifts of the data between -10:-5 sec and 5:10 sec
shift = round([-10:0.001:-5,5:0.001:10]*Fs);
shift_perm = randperm(length(shift),1000);
%Calculate the maximum cross correlation value of the shifted versions
Max_XCorr_shifted = zeros(1,1000);
for i=1:1000
    pad = zeros(abs(shift(shift_perm(i))),1);
    if shift(shift_perm(i))<0
        InstAmp2_shifted = [InstAmp2(abs(shift(shift_perm(i)))+1:end);pad];
    else
        InstAmp2_shifted = [pad;InstAmp2(1:length(InstAmp2)-shift(shift_perm(i)))];
    end
    AmpXCorr_shifted = xcorr(InstAmp1,InstAmp2_shifted,MaxLag,'coeff');
    Max_XCorr_shifted(i) = max(AmpXCorr_shifted);
end

% Determine the 95% confidence interval (mean + 2 standard deviations)
CI_95 = mean(Max_XCorr_shifted)+2*std(Max_XCorr_shifted);
% Result is considered significant if greater than 95% confidence interval
if max(AmpXCorr)>CI_95
    AmpXCorrSig = 1;
else
    AmpXCorrSig = 0;
end