%% 1. MAIN-FILE UNSUPERVISED ECG-BASED ATRIAL FIBRILLATION DETECTION
clear; close all; clc;
% dwtmode('zpd','nodisp'); % zero-padding for wavelet decomposition

% load the data: contains single-lead ECG data, the sampling frequency, and
% the clinical annotations (1 = normal, 2 = AF)
load ecg;

%% annotations

% Assume annotations is your 532500x1 data
% annotations = ...  <- Here you should insert the code that loads your data

% Sample frequency
fs = 250; %Hz

% Calculate the time vector in seconds
time = (0:length(annotations)-1)/fs;

% Create a new figure
figure;

% Plot the data
plot(time, annotations, 'LineWidth', 1.5);

% Use yticks to set the y-axis ticks to only be 1 and 2 
yticks([1 2])

% Name the ticks for better understanding
yticklabels({'Normal', 'Atrial Fibrillation'})

% Provide labels and title for the plot
xlabel('Time (s)')
ylabel('Heart Condition')
title('Heart Condition over Time')

% Optionally, if you want a grid
grid on;



%% 1.1 Preprocessing

slen = length(ecg);
t=[1:slen]/fs;
f = [0:slen/2]*fs/slen;
X = fft(ecg);
X1 = abs(X(1:slen/2+1));


%% 1.1.2 Artifact removal




% Derivative-based filter
% Choose your pole value
pole = -0.99;
a = [1, pole];
b = [1,-1];
ecg_filtered = filter(b, a, ecg);

% Notch filter
% Specify filter parameters
notch_freq = 50;  % Notch frequency
bw = 1;  % Bandwidth

% Design Notch filter
[b,a] = iirnotch(notch_freq/(fs/2), bw/(fs/2));

% Apply Notch filter
ecg_filtered = filter(b,a,ecg_filtered);

% Savitzky-Golay filter
order = 7;  % Order of polynomial fit
framelen = 15;  % Frame length

% Apply Savitzky-Golay filter
ecg_filtered = sgolayfilt(ecg_filtered, order, framelen);

% Moving average filter
windowSize = 2; % Change this to control the amount of smoothing
b = ones(1, windowSize)/windowSize;
ecg_filtered = conv(ecg_filtered, b, 'same');








% Define the interval with large artifacts 1
start_time_artifact = 556;  % start time of the artifact interval (in seconds)
end_time_artifact = 560.5;  % end time of the artifact interval (in seconds)

% Define the corresponding section of the original signal
noisy_ecg = ecg(fs*start_time_artifact+1:fs*end_time_artifact);

% Apply a bandpass filter to the artifact interval
low_cutoff = 0.2;  % Low cutoff frequency (in Hz)
high_cutoff = 70;  % High cutoff frequency (in Hz)
[b, a] = butter(3, [low_cutoff high_cutoff]/(fs/2), 'bandpass');
ecg_filtered_artifact = filter(b, a, noisy_ecg);

% Apply a notch filter to remove 50 Hz noise
notch_freq = 50;  % Notch frequency (in Hz)
bw = 1;  % Bandwidth
[b, a] = iirnotch(notch_freq/(fs/2), bw/(fs/2));
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Apply a low-pass filter to remove high frequency noise
high_cutoff = 20;  % High cutoff frequency (in Hz)
[b, a] = butter(3, high_cutoff/(fs/2), 'low');
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Derivative-based filter
% Choose your pole value
pole = -0.97;
a = [1, pole];
b = [1,-1];
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Concatenate the filtered artifact interval with the rest of the cleaned signal
ecg_filtered = [ecg_filtered(1:fs*start_time_artifact); ecg_filtered_artifact; ecg_filtered(fs*end_time_artifact+1:end)];




% Define the interval with large artifacts 2
start_time_artifact = 615;  % start time of the artifact interval (in seconds)
end_time_artifact = 620;  % end time of the artifact interval (in seconds)

% Define the corresponding section of the original signal
noisy_ecg = ecg(fs*start_time_artifact+1:fs*end_time_artifact);

% Apply a bandpass filter to the artifact interval
low_cutoff = 0.5;  % Low cutoff frequency (in Hz)
high_cutoff = 40;  % High cutoff frequency (in Hz)
[b, a] = butter(3, [low_cutoff high_cutoff]/(fs/2), 'bandpass');
ecg_filtered_artifact = filter(b, a, noisy_ecg);

% Apply a notch filter to remove 50 Hz noise
notch_freq = 50;  % Notch frequency (in Hz)
bw = 1;  % Bandwidth
[b, a] = iirnotch(notch_freq/(fs/2), bw/(fs/2));
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Apply a high-pass filter to remove baseline drift
high_cutoff = 0.5;  % High cutoff frequency (in Hz)
[b, a] = butter(3, high_cutoff/(fs/2), 'high');
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Apply a moving average filter to smooth the signal
windowSize = 10; % window size for the moving average filter
b = ones(1, windowSize) / windowSize;
ecg_filtered_artifact = filter(b, 1, ecg_filtered_artifact);


% Concatenate the filtered artifact interval with the rest of the cleaned signal
ecg_filtered = [ecg_filtered(1:fs*start_time_artifact); ecg_filtered_artifact; ecg_filtered(fs*end_time_artifact+1:end)];


% Define the interval with large artifacts 3
start_time_artifact = 634;  % start time of the artifact interval (in seconds)
end_time_artifact = 640;  % end time of the artifact interval (in seconds)

% Define the corresponding section of the original signal
noisy_ecg = ecg(fs*start_time_artifact+1:fs*end_time_artifact);

% Apply a bandpass filter to the artifact interval
low_cutoff = 0.2;  % Low cutoff frequency (in Hz)
high_cutoff = 40;  % High cutoff frequency (in Hz)
[b, a] = butter(3, [low_cutoff high_cutoff]/(fs/2), 'bandpass');
ecg_filtered_artifact = filter(b, a, noisy_ecg);

% Apply a notch filter to remove 50 Hz noise
notch_freq = 50;  % Notch frequency (in Hz)
bw = 1;  % Bandwidth
[b, a] = iirnotch(notch_freq/(fs/2), bw/(fs/2));
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Apply a high-pass filter to remove baseline drift
high_cutoff = 0.5;  % High cutoff frequency (in Hz)
[b, a] = butter(3, high_cutoff/(fs/2), 'high');
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Apply a low-pass filter to remove high frequency noise
high_cutoff = 40;  % High cutoff frequency (in Hz)
[b, a] = butter(3, high_cutoff/(fs/2), 'low');
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Derivative-based filter
% Choose your pole value
pole = -0.97;
a = [1, pole];
b = [1,-1];
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Concatenate the filtered artifact interval with the rest of the cleaned signal
ecg_filtered = [ecg_filtered(1:fs*start_time_artifact); ecg_filtered_artifact; ecg_filtered(fs*end_time_artifact+1:end)];


% Define the interval with large artifacts 3
start_time_artifact = 640;  % start time of the artifact interval (in seconds)
end_time_artifact = 643;  % end time of the artifact interval (in seconds)

% Define the corresponding section of the original signal
noisy_ecg = ecg(fs*start_time_artifact+1:fs*end_time_artifact);

% Apply a bandpass filter to the artifact interval
low_cutoff = 2;  % Low cutoff frequency (in Hz)
high_cutoff = 40;  % High cutoff frequency (in Hz)
[b, a] = butter(3, [low_cutoff high_cutoff]/(fs/2), 'bandpass');
ecg_filtered_artifact = filter(b, a, noisy_ecg);

% Apply a notch filter to remove 50 Hz noise
notch_freq = 50;  % Notch frequency (in Hz)
bw = 1;  % Bandwidth
[b, a] = iirnotch(notch_freq/(fs/2), bw/(fs/2));
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Apply a high-pass filter to remove baseline drift
high_cutoff = 0.5;  % High cutoff frequency (in Hz)
[b, a] = butter(3, high_cutoff/(fs/2), 'high');
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Derivative-based filter
% Choose your pole value
pole = -0.99;
a = [1, pole];
b = [1,-1];
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Invert the signal
ecg_filtered_artifact = -ecg_filtered_artifact;

% Concatenate the filtered artifact interval with the rest of the cleaned signal
ecg_filtered = [ecg_filtered(1:fs*start_time_artifact); ecg_filtered_artifact; ecg_filtered(fs*end_time_artifact+1:end)];

% Define the interval with large artifacts 4
start_time_artifact = 647;  % start time of the artifact interval (in seconds)
end_time_artifact = 650;  % end time of the artifact interval (in seconds)

% Define the corresponding section of the original signal
noisy_ecg = ecg(fs*start_time_artifact+1:fs*end_time_artifact);

% Apply a bandpass filter to the artifact interval
low_cutoff = 2;  % Low cutoff frequency (in Hz)
high_cutoff = 40;  % High cutoff frequency (in Hz)
[b, a] = butter(3, [low_cutoff high_cutoff]/(fs/2), 'bandpass');
ecg_filtered_artifact = filter(b, a, noisy_ecg);

% Apply a notch filter to remove 50 Hz noise
notch_freq = 50;  % Notch frequency (in Hz)
bw = 1;  % Bandwidth
[b, a] = iirnotch(notch_freq/(fs/2), bw/(fs/2));
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Apply a high-pass filter to remove baseline drift
high_cutoff = 0.5;  % High cutoff frequency (in Hz)
[b, a] = butter(3, high_cutoff/(fs/2), 'high');
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Derivative-based filter
% Choose your pole value
pole = -0.99;
a = [1, pole];
b = [1,-1];
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Invert the signal
ecg_filtered_artifact = -ecg_filtered_artifact;


% Concatenate the filtered artifact interval with the rest of the cleaned signal
ecg_filtered = [ecg_filtered(1:fs*start_time_artifact); ecg_filtered_artifact; ecg_filtered(fs*end_time_artifact+1:end)];

% Define the interval with large artifacts 5
start_time_artifact = 650;  % start time of the artifact interval (in seconds)
end_time_artifact = 652;  % end time of the artifact interval (in seconds)

% Define the corresponding section of the original signal
noisy_ecg = ecg(fs*start_time_artifact+1:fs*end_time_artifact);

% Apply a bandpass filter to the artifact interval
low_cutoff = 2;  % Low cutoff frequency (in Hz)
high_cutoff = 40;  % High cutoff frequency (in Hz)
[b, a] = butter(3, [low_cutoff high_cutoff]/(fs/2), 'bandpass');
ecg_filtered_artifact = filter(b, a, noisy_ecg);

% Apply a notch filter to remove 50 Hz noise
notch_freq = 50;  % Notch frequency (in Hz)
bw = 1;  % Bandwidth
[b, a] = iirnotch(notch_freq/(fs/2), bw/(fs/2));
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Apply a high-pass filter to remove baseline drift
high_cutoff = 0.5;  % High cutoff frequency (in Hz)
[b, a] = butter(3, high_cutoff/(fs/2), 'high');
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Derivative-based filter
% Choose your pole value
pole = -0.99;
a = [1, pole];
b = [1,-1];
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Invert the signal
ecg_filtered_artifact = -ecg_filtered_artifact;


% Concatenate the filtered artifact interval with the rest of the cleaned signal
ecg_filtered = [ecg_filtered(1:fs*start_time_artifact); ecg_filtered_artifact; ecg_filtered(fs*end_time_artifact+1:end)];


% Define the interval with large artifacts 6
start_time_artifact = 654,5;  % start time of the artifact interval (in seconds)
end_time_artifact = 656;  % end time of the artifact interval (in seconds)

% Define the corresponding section of the original signal
noisy_ecg = ecg(fs*start_time_artifact+1:fs*end_time_artifact);

% Apply a bandpass filter to the artifact interval
low_cutoff = 2;  % Low cutoff frequency (in Hz)
high_cutoff = 40;  % High cutoff frequency (in Hz)
[b, a] = butter(3, [low_cutoff high_cutoff]/(fs/2), 'bandpass');
ecg_filtered_artifact = filter(b, a, noisy_ecg);

% Apply a notch filter to remove 50 Hz noise
notch_freq = 50;  % Notch frequency (in Hz)
bw = 1;  % Bandwidth
[b, a] = iirnotch(notch_freq/(fs/2), bw/(fs/2));
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Apply a high-pass filter to remove baseline drift
high_cutoff = 0.5;  % High cutoff frequency (in Hz)
[b, a] = butter(3, high_cutoff/(fs/2), 'high');
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Derivative-based filter
% Choose your pole value
pole = -0.99;
a = [1, pole];
b = [1,-1];
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Invert the signal
ecg_filtered_artifact = -ecg_filtered_artifact;


% Concatenate the filtered artifact interval with the rest of the cleaned signal
ecg_filtered = [ecg_filtered(1:fs*start_time_artifact); ecg_filtered_artifact; ecg_filtered(fs*end_time_artifact+1:end)];


% Define the interval with large artifacts 7
start_time_artifact = 657.5;  % start time of the artifact interval (in seconds)
end_time_artifact = 659;  % end time of the artifact interval (in seconds)

% Define the corresponding section of the original signal
noisy_ecg = ecg(fs*start_time_artifact+1:fs*end_time_artifact);

% Apply a bandpass filter to the artifact interval
low_cutoff = 0.5;  % Low cutoff frequency (in Hz)
high_cutoff = 40;  % High cutoff frequency (in Hz)
[b, a] = butter(3, [low_cutoff high_cutoff]/(fs/2), 'bandpass');
ecg_filtered_artifact = filter(b, a, noisy_ecg);

% Apply a notch filter to remove 50 Hz noise
notch_freq = 50;  % Notch frequency (in Hz)
bw = 1;  % Bandwidth
[b, a] = iirnotch(notch_freq/(fs/2), bw/(fs/2));
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Apply a high-pass filter to remove baseline drift
high_cutoff = 2;  % High cutoff frequency (in Hz)
[b, a] = butter(3, high_cutoff/(fs/2), 'high');
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Apply a moving average filter to smooth the signal
windowSize = 10; % window size for the moving average filter
b = ones(1, windowSize) / windowSize;
ecg_filtered_artifact = filter(b, 1, ecg_filtered_artifact);

% Derivative-based filter
% Choose your pole value
pole = -0.99;
a = [1, pole];
b = [1,-1];
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Invert the signal
ecg_filtered_artifact = -ecg_filtered_artifact;


% Concatenate the filtered artifact interval with the rest of the cleaned signal
ecg_filtered = [ecg_filtered(1:fs*start_time_artifact); ecg_filtered_artifact; ecg_filtered(fs*end_time_artifact+1:end)];


% Define the interval with large artifacts 8
start_time_artifact = 674;  % start time of the artifact interval (in seconds)
end_time_artifact = 676;  % end time of the artifact interval (in seconds)

% Define the corresponding section of the original signal
noisy_ecg = ecg(fs*start_time_artifact+1:fs*end_time_artifact);

% Apply a bandpass filter to the artifact interval
low_cutoff = 2;  % Low cutoff frequency (in Hz)
high_cutoff = 50;  % High cutoff frequency (in Hz)
[b, a] = butter(3, [low_cutoff high_cutoff]/(fs/2), 'bandpass');
ecg_filtered_artifact = filter(b, a, noisy_ecg);

% Apply a moving average filter to smooth the signal
windowSize = 10; % window size for the moving average filter
b = ones(1, windowSize) / windowSize;
ecg_filtered_artifact = filter(b, 1, ecg_filtered_artifact);

% Derivative-based filter
% Choose your pole value
pole = -0.93;
a = [1, pole];
b = [1,-1];
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Invert the signal
ecg_filtered_artifact = -ecg_filtered_artifact;

% Concatenate the filtered artifact interval with the rest of the cleaned signal
ecg_filtered = [ecg_filtered(1:fs*start_time_artifact); ecg_filtered_artifact; ecg_filtered(fs*end_time_artifact+1:end)];


% Define the interval with large artifacts 8
start_time_artifact = 738.5;  % start time of the artifact interval (in seconds)
end_time_artifact = 739;  % end time of the artifact interval (in seconds)

% Define the corresponding section of the original signal
noisy_ecg = ecg(fs*start_time_artifact+1:fs*end_time_artifact);

% Apply a bandpass filter to the artifact interval
low_cutoff = 2;  % Low cutoff frequency (in Hz)
high_cutoff = 50;  % High cutoff frequency (in Hz)
[b, a] = butter(3, [low_cutoff high_cutoff]/(fs/2), 'bandpass');
ecg_filtered_artifact = filter(b, a, noisy_ecg);

% Apply a moving average filter to smooth the signal
windowSize = 10; % window size for the moving average filter
b = ones(1, windowSize) / windowSize;
ecg_filtered_artifact = filter(b, 1, ecg_filtered_artifact);

% Derivative-based filter
% Choose your pole value
pole = -0.93;
a = [1, pole];
b = [1,-1];
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Invert the signal
ecg_filtered_artifact = -ecg_filtered_artifact;

% Concatenate the filtered artifact interval with the rest of the cleaned signal
ecg_filtered = [ecg_filtered(1:fs*start_time_artifact); ecg_filtered_artifact; ecg_filtered(fs*end_time_artifact+1:end)];

% Define the interval with large artifacts 8
start_time_artifact = 741;  % start time of the artifact interval (in seconds)
end_time_artifact = 743;  % end time of the artifact interval (in seconds)

% Define the corresponding section of the original signal
noisy_ecg = ecg(fs*start_time_artifact+1:fs*end_time_artifact);

% Apply a bandpass filter to the artifact interval
low_cutoff = 1;  % Low cutoff frequency (in Hz)
high_cutoff = 50;  % High cutoff frequency (in Hz)
[b, a] = butter(3, [low_cutoff high_cutoff]/(fs/2), 'bandpass');
ecg_filtered_artifact = filter(b, a, noisy_ecg);

% Apply a moving average filter to smooth the signal
windowSize = 5; % window size for the moving average filter
b = ones(1, windowSize) / windowSize;
ecg_filtered_artifact = filter(b, 1, ecg_filtered_artifact);

% Derivative-based filter
% Choose your pole value
pole = -0.93;
a = [1, pole];
b = [1,-1];
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Invert the signal
ecg_filtered_artifact = -ecg_filtered_artifact;

% Concatenate the filtered artifact interval with the rest of the cleaned signal
ecg_filtered = [ecg_filtered(1:fs*start_time_artifact); ecg_filtered_artifact; ecg_filtered(fs*end_time_artifact+1:end)];


% Define the interval with large artifacts 8
start_time_artifact = 748;  % start time of the artifact interval (in seconds)
end_time_artifact = 750;  % end time of the artifact interval (in seconds)

% Define the corresponding section of the original signal
noisy_ecg = ecg(fs*start_time_artifact+1:fs*end_time_artifact);

% Apply a bandpass filter to the artifact interval
low_cutoff = 1;  % Low cutoff frequency (in Hz)
high_cutoff = 50;  % High cutoff frequency (in Hz)
[b, a] = butter(3, [low_cutoff high_cutoff]/(fs/2), 'bandpass');
ecg_filtered_artifact = filter(b, a, noisy_ecg);

% Apply a moving average filter to smooth the signal
windowSize = 5; % window size for the moving average filter
b = ones(1, windowSize) / windowSize;
ecg_filtered_artifact = filter(b, 1, ecg_filtered_artifact);

% Derivative-based filter
% Choose your pole value
pole = -0.93;
a = [1, pole];
b = [1,-1];
ecg_filtered_artifact = filter(b, a, ecg_filtered_artifact);

% Invert the signal
ecg_filtered_artifact = -ecg_filtered_artifact;

% Concatenate the filtered artifact interval with the rest of the cleaned signal
ecg_filtered = [ecg_filtered(1:fs*start_time_artifact); ecg_filtered_artifact; ecg_filtered(fs*end_time_artifact+1:end)];





% plot 1: original ecg and power spectrum of original ecg + filtered sig. 

% plot the original ecg + filtered signal
figure
subplot(2,1,1)
plot(t,ecg)
hold on
plot(t,ecg_filtered)
xlabel('Seconds')
ylabel('Signal')
title('Filtered signal: Complete signal')
legend('Original ECG','Filtered ECG') 
% Set the x-axis and y-axis limits of the plot
xlim([0 2130]);
ylim([-4 4]);


y1 = fft(ecg_filtered);


subplot(2,1,2)
plot(f,20*log10(X1))
hold on
plot(f,20*log10(abs(y1(1:slen/2+1))))
xlabel('Hertz')
ylabel('Magnitude spectrum')
title('Spectrum of the filtered signal')
legend('Original ECG','Filtered ECG')

% plot 2: original ecg and power spectrum of original ecg + filtered sig. 

% plot the original ecg + filtered signal
figure

plot(t,ecg)
hold on
plot(t,ecg_filtered)
xlabel('Seconds')
ylabel('Signal')
title('Filtered signal: Large artefact')
legend('Original ECG','Filtered ECG') 
% Set the x-axis and y-axis limits of the plot
xlim([633 652]);
ylim([-5.5 5]);


% plot 3: original ecg and power spectrum of original ecg + filtered sig. 

% plot the original ecg + filtered signal
figure

plot(t,ecg)
hold on
plot(t,ecg_filtered)
xlabel('Seconds')
ylabel('Signal')
title('Filtered signal: Large artefact')
legend('Original ECG','Filtered ECG') 
% Set the x-axis and y-axis limits of the plot
xlim([612 624]);
ylim([-2.5 5]);


% plot 4: original ecg and power spectrum of original ecg + filtered sig. 

% plot the original ecg + filtered signal
figure

plot(t,ecg)
hold on
plot(t,ecg_filtered)
xlabel('Seconds')
ylabel('Signal')
title('Filtered signal: Large artefact')
legend('Original ECG','Filtered ECG') 

% Set the x-axis and y-axis limits of the plot
xlim([548 569]);
ylim([-4 5]);

% plot 5: original ecg and power spectrum of original ecg + filtered sig. 

% plot the original ecg + filtered signal
figure

plot(t,ecg)
hold on
plot(t,ecg_filtered)
xlabel('Seconds')
ylabel('Signal')
title('Filtered signal: Large artefact')
legend('Original ECG','Filtered ECG') 

% Set the x-axis and y-axis limits of the plot
xlim([673 676]);
ylim([-1 4]);


% plot 6: original ecg and power spectrum of original ecg + filtered sig. 

% plot the original ecg + filtered signal
figure

plot(t,ecg)
hold on
plot(t,ecg_filtered)
xlabel('Seconds')
ylabel('Signal')
title('Filtered signal: Large artefact')
legend('Original ECG','Filtered ECG') 

% Set the x-axis and y-axis limits of the plot
xlim([738 743]);
ylim([-2 5]);

% plot 7: original ecg and power spectrum of original ecg + filtered sig. 

% plot the original ecg + filtered signal
figure

plot(t,ecg)
hold on
plot(t,ecg_filtered)
xlabel('Seconds')
ylabel('Signal')
title('Filtered signal: Large artefact')
legend('Original ECG','Filtered ECG') 

% Set the x-axis and y-axis limits of the plot
xlim([747 751]);
ylim([-2 4.5]);

% plot 8: original ecg and power spectrum of original ecg + filtered sig. 

% plot the original ecg + filtered signal
figure

plot(t,ecg)
hold on
plot(t,ecg_filtered)
xlabel('Seconds')
ylabel('Signal')
title('Filtered signal: Normal segment')
legend('Original ECG','Filtered ECG') 

% Set the x-axis and y-axis limits of the plot
xlim([1500 1510]);
ylim([-2 4.5]);

% plot 9: original ecg and power spectrum of original ecg + filtered sig. 

% plot the original ecg + filtered signal
figure

plot(t,ecg)
hold on
plot(t,ecg_filtered)
xlabel('Seconds')
ylabel('Signal')
title('Filtered signal: Normal segment')
legend('Original ECG','Filtered ECG') 

% Set the x-axis and y-axis limits of the plot
xlim([1500 1510]);
ylim([-2 4.5]);
%% 1.2 Segmentation

% Use the buffer function to segment the signal into 10-second segments with no overlap


% check: do you have one label per 10s-segment?

% Segment the ECG signal into non-overlapping segments of 10 seconds
segments = buffer(ecg_filtered, 2500);

% Divide annotations vector into 213 segments
segments_annotations = reshape(annotations, [], 213);



figure('Name','Segmentation');
subplot(2,1,1);

slen_segment_normal=length(segments(:,1))
t_segment=[1:slen_segment_normal]/fs+70;

% Normal segment

plot(t_segment,ecg_filtered(17501:20000))
xlabel('Seconds')
ylabel('Signal')
title('Normal segment')

slen_segment_AF=length(segments(:,1))
t_segment=[1:slen_segment_AF]/fs+20;

slen_segment_AF=length(segments(:,1))
t_segment=[1:slen_segment_AF]/fs+20;

subplot(2,1,2);
% AF segment
plot(t_segment,ecg_filtered(5001:7500))
xlabel('Seconds')
ylabel('Signal')
title('AF segment')

%% 1.3 Feature extraction

% check: do you have a N x 16 or 16 x N feature matrix, with N the number
% of segments? yes 213x16 matrix

% The calculation of the features are numbered 1-6

% Calculate the signal energies of each segment

% Calculate the squared matrix
squared_matrix = segments.^2;

% Calculate the mean-squared value of each segment
signal_energies = sum(squared_matrix, 1);






% 2) Calculate the wavelet decomposition of the ECG signal using Daubechies


% Define the wavelet to use for the decomposition
wavelet = 'db4';

% Define the level of decomposition to use
level = 5;

% Initialize a matrix to store the wavelet energies of each segment
waveletEnergies = zeros(level+1, size(segments, 2));

% Initialize a matrix to store the relative wavelet energies of each segment
relativeWaveletEnergies = zeros(level+1, size(segments, 2));

% Iterate through the columns of the matrix
for i = 1:size(segments, 2)
    % Perform the wavelet decomposition on the current segment
    [c, l] = wavedec(segments(:, i), level, wavelet);
    
    % Compute the energies for the detail coefficients at each level
    for j = 1:level
        % Extract the detail coefficients at level j
        cD = detcoef(c, l, j);
        
        % Compute the energy at this level
        waveletEnergies(j, i) = sum(cD.^2);
    end
    
    % Extract the final approximation coefficients
    cA = appcoef(c, l, wavelet, level);
    
    % Compute the energy of the approximation coefficients
    waveletEnergies(level+1, i) = sum(cA.^2);

    % Compute the total energy
    totalEnergy = sum(waveletEnergies(:, i));

    % Compute the relative energies
    relativeWaveletEnergies(:, i) = waveletEnergies(:, i) / totalEnergy;
end












% Finding the R peaks and visualizing on some segments
% Create an empty cell array to store the R peak amplitudes for each segment
amplitudes = {};

for i = 1:213
    % Use the findpeaks function to detect the R peaks in the segment
    [pks, locs] = findpeaks(segments(:,i), 'MinPeakHeight', 1.5*std(segments(:,i)));

    % Convert the peak locations to seconds
    segmentTimeInSeconds = locs / fs;

    % Add the peak locations and amplitudes from the current segment to the appropriate cell arrays
    timeInSeconds{i} = segmentTimeInSeconds;
    amplitudes{i} = pks;


end

% Create an empty vector to store the number of R peaks in each segment
number_peaks_segment = zeros(1, length(timeInSeconds));

% Loop through each segment
for i = 1:length(timeInSeconds)
    % Get the number of elements in the current vector in the timeInSeconds cell array
    numElements = length(timeInSeconds{i});

    % Store the number of elements in the number_peaks_segment vector
    number_peaks_segment(i) = numElements;
end


% Create an empty vector to store the number of RR intervals in each segment
number_rr_segment = zeros(1, length(timeInSeconds));

% Loop through each segment
for i = 1:length(timeInSeconds)
    % Get the vector of R peak locations for the current segment
    segmentTimeInSeconds = timeInSeconds{i};

    % Use the diff function to compute the differences between consecutive elements in the vector of R peak locations
    rrIntervals = diff(segmentTimeInSeconds);

    % The number of RR intervals in the segment is equal to the number of elements in the rrIntervals vector, minus one
    numRRIntervals = length(rrIntervals) - 1;

    % Store the number of RR intervals in the number_rr_segment vector
    number_rr_segment(i) = numRRIntervals;
end


% Visualizing the R peak detection on some segments: 

% Normal segment
subplot(2,1,1);
plot(t_segment,ecg_filtered(1:2500))
hold on

% Loop through the columns of timeInSeconds
for i = 1:size(timeInSeconds, 2)
    % Plot vertical lines at the time of each R peak
    for j = 1:length(timeInSeconds{i})
        plot([timeInSeconds{i}(j) timeInSeconds{i}(j)], ylim, 'r')
    end
end

xlabel('Seconds')
ylabel('Signal')
title('Normal segment')



% 4) Create an empty vector to store the RMSSD for each segment
rmssd_segment = zeros(1, length(timeInSeconds));

% Loop through each segment
for i = 1:length(timeInSeconds)
    % Get the vector of R peak locations for the current segment
    segmentTimeInSeconds = timeInSeconds{i};

    % Use the diff function to compute the differences between consecutive elements in the vector of R peak locations
    rrIntervals = diff(segmentTimeInSeconds);

    % Square each element in the rrIntervals vector
    squaredDifferences = rrIntervals.^2;

    % Take the mean of the squared differences
    meanSquaredDifferences = mean(squaredDifferences);

    % Take the square root of the mean of the squared differences to get the RMSSD for the segment
    rmssd = sqrt(meanSquaredDifferences);

    % Store the RMSSD for the current segment in the rmssd_segment vector
    rmssd_segment(i) = rmssd;
end


% Initialize an empty cell array to hold the results
fractions = cell(1, 213);





% 5) Fraction of successive RR intervals that differ by more than 50ms for each segment
% Loop over each segment in timeInSeconds
for i = 1:213
    
    % Compute the differences between successive elements in the current segment
    rrIntervals = diff(timeInSeconds{i});
    
    % Find the indices of the elements in rrIntervals that have a magnitude greater than 0.05
    indices = find(abs(rrIntervals) > 0.05);
    
    % Calculate the fraction of elements with a magnitude greater than 0.05
    fractions{i} = numel(indices) / numel(timeInSeconds{i});
end

% Convert the fractions cell array into a row vector
fractionsVector = cell2mat(fractions);





% 6) Calculating the minimal RR-interval (minRR)

rrIntervals = cellfun(@diff, timeInSeconds, 'UniformOutput', false);

% Apply the min function to each element of the resulting cell array
minRR = cellfun(@min, rrIntervals, 'UniformOutput', false);

if isempty(minRR{57})
    minRR{57} = 0;
end

dim = size(minRR);
minRRRowvector = zeros(dim);

% Convert the cell array into a matrix
minRRRowvector = cell2mat(minRR);




% making a 16x213 matrix from all the features

featurematrix=[signal_energies;waveletEnergies;relativeWaveletEnergies;rmssd_segment;fractionsVector;minRRRowvector]

% Identify the missing or invalid values in the 57th column
missing = isnan(featurematrix(:, 57)) | isinf(featurematrix(:, 57));

% Replace the missing or invalid values with zeros
featurematrix(missing, 57) = 0;

std_segments = std(featurematrix,0,2);
mean_segments = mean(featurematrix,2);

% Create the table
table_mean_std = table(mean_segments, std_segments);

% Transposing feature matrix ->  making ready for dimensionality reduction
featurematrix = featurematrix.'


figure('Name','Feature extraction')


%% Tachogram + RR Interval Visualization

figure

subplot(2,2,1)
% Normal R-peak detection

% 8th normal segment

slen_segment_normal=length(segments(:,1))
t_segment=[1:slen_segment_normal]/fs+70;

plot(t_segment, ecg_filtered(17501:20000))
hold on;

% Use the 7th cell of the timeInSeconds array for the x values of the red line
rPeaks = timeInSeconds{8}+70;

% Find the y values of the red line at the locations specified by the timeInSeconds vector
ecgPeaks = ecg_filtered(17500+timeInSeconds{8} * fs);

scatter(rPeaks, ecgPeaks, 25, 'r');
xlabel('Seconds')
ylabel('Signal')
title('8th segment (normal)')



subplot(2,2,2)
% Normal tachogram

% Convert R-peak locations from seconds to samples
rPeaks8 = cellfun(@(x) x * fs, timeInSeconds,'UniformOutput',false);

% Get R-peak locations for segment 3
rPeaksSegment8 = rPeaks8{8};

% Calculate tachogram values for segment 3
tachogramSegment8 = diff(rPeaksSegment8) * 1000 / fs;

% Plot tachogram for segment 3
plot(linspace(70, 80, length(tachogramSegment8)), tachogramSegment8)

% Set x-axis labels to show time in seconds
xlabel('Time (s)')

% Set y-axis labels to show tachogram values
ylabel('RR interval (ms)')

% Set x-axis labels to show time in seconds
xlabel('Time (s)')
title('Tachogram 8th segment (normal)')






subplot(2,2,3)
% AF R-peak detection
slen_segment_normal=length(segments(:,1))
t_segment=[1:slen_segment_normal]/fs+20;

plot(t_segment, ecg_filtered(5001:7500))
hold on;

% Use the 8th cell of the timeInSeconds array for the x values of the red line
rPeaks = timeInSeconds{3}+20;

% Find the y values of the red line at the locations specified by the timeInSeconds vector
ecgPeaks = ecg_filtered(5000+timeInSeconds{3} * fs);

scatter(rPeaks, ecgPeaks, 25, 'r');
xlabel('Seconds')
ylabel('Signal')
title('3th segment (AF)')


subplot(2,2,4)
% AF tachogram

% Convert R-peak locations from seconds to samples
rPeaks3 = cellfun(@(x) x * fs, timeInSeconds,'UniformOutput',false);

% Get R-peak locations for segment 3
rPeaksSegment3 = rPeaks3{3};

% Calculate tachogram values for segment 3
tachogramSegment3 = diff(rPeaksSegment3) * 1000 / fs;

% Plot tachogram for segment 3
plot(linspace(20, 30, length(tachogramSegment3)), tachogramSegment3)

% Set x-axis labels to show time in seconds
xlabel('Time (s)')

% Set y-axis labels to show tachogram values
ylabel('RR interval (ms)')

% Set x-axis labels to show time in seconds
xlabel('Time (s)')
title('Tachogram 3th segment (AF)')



% Visualizing the R peak detection on some segments:

% 3th AF segment
figure('Name','Segmentation');
subplot(5,1,1);

slen_segment_normal=length(segments(:,1))
t_segment=[1:slen_segment_normal]/fs+20;

plot(t_segment, ecg_filtered(5001:7500))
hold on;

% Use the 3rd cell of the timeInSeconds array for the x values of the red line
rPeaks = timeInSeconds{3}+20;

% Find the y values of the red line at the locations specified by the timeInSeconds vector
ecgPeaks = ecg_filtered(5000+timeInSeconds{3} * fs);


scatter(rPeaks, ecgPeaks, 25, 'r');
xlabel('Seconds')
ylabel('Signal')
title('3th AF segment')





% 4th AF segment
subplot(5,1,2);

slen_segment_normal=length(segments(:,1))
t_segment=[1:slen_segment_normal]/fs+30;


plot(t_segment, ecg_filtered(7501:10000))
hold on;

% Use the 8th cell of the timeInSeconds array for the x values of the red line
rPeaks = timeInSeconds{4}+30;

% Find the y values of the red line at the locations specified by the timeInSeconds vector
ecgPeaks = ecg_filtered(7500+timeInSeconds{4} * fs);

scatter(rPeaks, ecgPeaks, 25, 'r');
xlabel('Seconds')
ylabel('Signal')
title('4th AF segment')






% 8th normal segment

subplot(5,1,4);

slen_segment_normal=length(segments(:,1))
t_segment=[1:slen_segment_normal]/fs+70;

plot(t_segment, ecg_filtered(17501:20000))
hold on;

% Use the 7th cell of the timeInSeconds array for the x values of the red line
rPeaks = timeInSeconds{8}+70;

% Find the y values of the red line at the locations specified by the timeInSeconds vector
ecgPeaks = ecg_filtered(17500+timeInSeconds{8} * fs);

scatter(rPeaks, ecgPeaks, 25, 'r');
xlabel('Seconds')
ylabel('Signal')
title('8th normal segment')



% 11th normal segment

subplot(5,1,5);

slen_segment_normal=length(segments(:,1))
t_segment=[1:slen_segment_normal]/fs+100;

plot(t_segment, ecg_filtered(25001:27500))
hold on;

% Use the 7th cell of the timeInSeconds array for the x values of the red line
rPeaks = timeInSeconds{11}+100;

% Find the y values of the red line at the locations specified by the timeInSeconds vector
ecgPeaks = ecg_filtered(25000+timeInSeconds{11} * fs);

scatter(rPeaks, ecgPeaks, 25, 'r');
xlabel('Seconds')
ylabel('Signal')
title('11th normal segment')




    
%% 1.4 Feature visualization and analysis

% 1.4.1 Normalization
% Units of your features are not the same --> ensure all features are given equal weight in the clustering process, improves the accuracy of the model
% first row --> large numbers --> might affect covariance matrix 


% Normalize the data in the feature matrix using zscore 
featurematrix_normalized = zscore(featurematrix);


% 1.4.2 Dimensionality reduction
% Dimensionality reduction -> reduces #features -> easier for classification

% Perform principal component analysis (PCA) on the standardized features
% score -> projections of orginal data points onto principal components ->
% represent data in lower-dimensional space
[coeff, score, latent] = pca(featurematrix_normalized);

% Retain only the first two principal components
coeff = coeff(:, 1:8);
score = score(:, 1:8);

% Transform the features using the retained principal components
% multiply original feature matrix by coeff matrix to create new matrix
% features_transformed, projecting original feature matrix onto space
% spanned by first data components = reduced representation of original
% data
features_transformed = featurematrix_normalized * coeff;

% The eigenvalues of the covariance

% Compute the eigenvalues of the normalized feature matrix
% use normalized version of the matrix -> ensures features on same scale
cov_matrix = cov(featurematrix_normalized);

% Compute the eigenvalues and eigenvectors of the covariance matrix
[eigenvectors, eigenvalues] = eig(cov_matrix);

% Plot the eigenvalues
% 16 eigenvalues corresponding to the 16 features of the featurematrix
figure('Name','Eigenvalues covariance')
plot(sort(diag(eigenvalues), 'descend'))

figure('Name','Eigenvalues covariance')
plot(sort(diag(eigenvalues), 'descend'))

% 1.4.3 Plot the reduced feature space

% Create a new figure for the plot
figure('Name','Reduced feature space')

% Plot the transformed features in the reduced two-dimensional space
% scatter plot will have 16 points, each representing a row in the feature matrix
% number of columns in the feature matrix does not affect the number of points in the scatter plot

scatter(score(:,1), score(:,2));

% Add labels to the x- and y-axes
xlabel('Feature 1');
ylabel('Feature 2');

% 1.4.4 feature matrix rank defficient?
% rank of feature matrix = 16 -> not rank defficient -> all rows are
% linearily independent

% Perform singular value decomposition on the matrix
[U,S,V] = svd(featurematrix);

% Check the number of non-zero singular values
num_nonzero_singular_values = sum(diag(S) > 0);

% Print the rank of the matrix
fprintf('The rank of the feature matrix is: %d\n', num_nonzero_singular_values);

%% 1.5 Unsupervised classification

% 1.5.1 Evaluation

% Calculate the pairwise Manhattan distances between the points using the pdist2 function
D = pdist2(score, score, 'cityblock');

% Perform K-means clustering
opts = statset('Display','final');
[clusters, centroids] = kmeans(score, 2, 'Distance','sqeuclidean', 'Start','sample', 'Replicates',5, 'Options',opts);

% Create a new figure for the plot
figure('Name','K-means clustering')

% Add the cluster centroids to the plot
hold on

colormap(cool(length(unique(clusters))))
scatter(score(:,1), score(:,2), [], clusters)

scatter(centroids(:,1), centroids(:,2), 'filled', 'MarkerFaceColor', 'k')

hold off


% accuracy, specificity and 

% Extract the first row of the matrix, which contains the first element from each column
labels = segments_annotations(1,:).';

% Compute the confusion matrix
confusion_matrix = confusionmat(labels, clusters);

% Extract the true positive, true negative, false positive, and false negative counts
tp = confusion_matrix(1,1); % true positives
tn = confusion_matrix(2,2); % true negatives
fp = confusion_matrix(1,2); % false positives
fn = confusion_matrix(2,1); % false negatives

% Compute the accuracy, sensitivity and specificity
accuracy = (tp + tn) / (tp + tn + fp + fn);
sensitivity = tp / (tp + fn);
specificity = tn / (tn + fp);

% Create a table containing the accuracy, sensitivity, and specificity
results_table = array2table([accuracy, sensitivity, specificity], 'VariableNames', {'Accuracy', 'Sensitivity', 'Specificity'});

% Define the vertices of the region of interest
x = [-30 20 -30 30];
y = [-10 -10 30 30];

% Use the inpolygon() function to find the points inside the region
inside = inpolygon(score(:,1), score(:,2), x, y);

% Count the number of points inside the region
num_points = sum(inside);
