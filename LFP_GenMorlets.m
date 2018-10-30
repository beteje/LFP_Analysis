function W = LFP_GenMorlets(F,fs)
%% Generate Morlet wavelets
% 
%   Inputs
%   F   -- vector of central frequencies [Hz]
%   fs  -- sampling frequency
%
%   Output
%   W   -- cell of vectors containing coefficients for each wavelet

W   =   cell(1,length(F)); 
ts  =   1/fs;   % sample period, second
nc  =   4;      % number of cycles in each wavelet

for i = 1:length(F)
    f0 = F(i); % select central frequency for new wavelet
    
    sigma_f = f0/nc; % standard deviation in frequency
    sigma_t = 1/(2*pi*sigma_f); % standard deviation in time
    
    t = 0:ts:5*sigma_t; % time vector over which to calculate wavelet coefficients
    t = [-t(end:-1:2) t];
    
    % generate wavelet coefficients for frequency f0
    W{i} = (sigma_t*sqrt(pi))^-0.5 * exp(-t.^2/(2*sigma_t^2)) .* exp(1i*2*pi*f0*t);
end
