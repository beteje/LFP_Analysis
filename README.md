# LFP_Analysis

All code can be run from the Start_LFP_Analysis.m file
All the parameters which need to be set for subsequent analysis are contained in this file

1.	First choose which analysis to perform – for each parameter set to 1 if you want to perform the analysis and 0 if you don’t want to perform the analysis

    •	PhaseLocking
    
    •	SpikeFieldCoherence
    
    •	OffsetPhaseLocking
    
    •	AmplitudeXCorr
    
2.	Next choose whether to perform single spike train to its corresponding LFP channel or for the LFP to LFP analysis one channel in one region to one channel in the other region. Or perform single spike train to multiple LFP channels or for the LFP to LFP analysis all channels in one region to all channels in the other region.

    •	LFPSpike – 1 = one to one, 0 = one to multiple
    
The following parameters are set according to the data you are analyzing and what settings you want.

3.	The first few are generic parameters no matter which analysis you are performing

    •	DataFile – the name of the file where your data is stored
    
    •	FileName – the name of the file where you want to save the results
    
    •	NumChannels  – the number of channels which have been recorded (the code assumes that you have the same number of channels in each region (i.e. NumChannels/2)
    
    •	LFPRate – the sampling rate of the LFP data in hertz
    
    •	FiltFreq – the cut off frequency of the low pass filter for the LFP data in hertz
    
    •	MinNumSpikes – the minimum number of spikes which we require for analysis (to ensure statistical significance of results
    
    •	AnalyRange – start and end times of the blocks of data you want to perform the analysis on in seconds
    
4.	Parameters for the phase locking analysis

    •	PL_WavFreq – the frequencies you want to perform the phase locking analysis at – these are used to create the wavelets in hertz
    
5.	Parameters for the spike field coherence analysis

    •	SFC_SegLength – the length of the segments either side of the spikes to be used for the SFC analysis in milliseconds. The total segment length will be 2*SFC_SegLength (i.e. the spike time ±SFC_SegLength).
    
6.	Parameters for the offset phase locking analysis (only need to be set if using LFPSpike=1 otherwise the analysis is not run)

    •	OPL_FreqBand – the minimum and maximum frequencies of the frequency band you want to analyze the phase locking over in hertz
    
    •	OPL_MaxLag – the maximum time offset to be analyzed in milliseconds
    
    •	OPL_StepSize – the time steps to analyze the phase locking at in milliseconds
    
7.	Parameters for the amplitude cross correlation analysis

    •	XCorr_MaxLag – maximum lag for cross correlation in milliseconds
    
    •	XCorr_FreqBands – maximum and minimum frequencies of the frequency bands to analyze the amplitude cross correlation within in hertz
    
    •	XCorr_AnalyRange – the start and end times of the blocks of data you want to perform the amplitude cross correlation on in seconds
    
    •	XCorr_Channels – pairs of channels you want to compare the amplitude cross correlation between (only needs to be set if using LFPSpike=1 otherwise all channels are used)
    
After setting all of the parameters run this code and it will store all the parameters in a structure and pass them to the analysis code. The code calls either LFP_Analysis.m or LFP_Analysis_XChannel.m depending on the setting of the LFPSpike parameter.
