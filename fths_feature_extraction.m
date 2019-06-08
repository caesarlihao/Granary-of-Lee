% This program extracts the 1/f feature angle of EEG signal.
% The specific algorithm is as follows:
% £¨1£©Load the .mat file and make the EEG data matrix in the form of M * N, where M is the number of electrodes, N is time-based EEG data.
% £¨2£©According to the segment length of 500ms, the EEG data is divided into non-overlapping segments,
%      and the last segment is filled with 0.
% £¨3£©Calculating power spectral density £¨PSD£© of each segment by FFT.
% £¨4£©The logarithms of PSD (y-axis) and frequency (x-axis) are calculated respectively,
%      and the 1/f relationship graph is obtained.
% £¨5£©Fitting the 1/f relation graph to get the 1/f feature straight line.
% £¨6£©The slope of the 1/f feature straight line is inversely triangulated to obtain 1/f feature angle.
% ¡ª¡ª Porgrammed by Li Hao, 2018.07


clear;
clc;

% Replace xxxx with your .mat or .dat file
% Suppose that the matrix storing EEG data is called eegdata
load xxxx.mat;

% Define sampling rate
% Replace ??? with your sampling rate in the EEG acquisition
fs=???;

% Define segment length (ms)
section_length=500;

% All EEG data are segmented according to section_length, and the last segment is filled with 0
section_num=ceil(length(eegdata)/fs*1000/section_length);
eegdata(:,length(eegdata)+1:section_num*fs*section_length/1000)=0;

% Define FFT points
fft_num=yyyy;

% Segments as the first level loops
for i=1:section_num
    
    % Electrodes as the second level loops
    for j=1:zzzz
        
        eeg_psd(j,:,i)=abs(fft(eegdata(j,(i-1)*fs*section_length/1000+1:i*fs*section_length/1000),fft_num));
        
        
        % Since the frequencies of EEG active components are usually 0.5-50 Hz, other frequency components are discarded.
        signal(j,:,i)=log(eeg_psd(j,ceil(0.5*fft_num/fs):ceil(50*fft_num/fs),i));
        
        % Calculate the logarithmic frequency as the x-axis
        f=log((ceil(0.5*fft_num/fs):(length(signal(j,:,i))+2))./(length(signal(j,:,i))+2).*50);
        
        % Calculate fitting parameter of a segment
        polyfit_para(j,:,i)=polyfit(f,signal(j,:,i),1);
        
        % Calculating the 1/f feature angle of the segment
        angle(j,i)=atan(polyfit_para(j,1,i))/pi*180;
    end;
end;