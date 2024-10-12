%% 2. MAIN FILE SUPERVISED EEG-BASED SLEEP STAGING
% Path initialisation
clear; close all; clc
addpath('./FastICA_25')
addpath('./Data');
%% Load the data
% sleep_staging.mat contains the following variables:
% channels: the channel names
% eeg: the eeg data, split in training data (eeg.train) and test data (eeg.test)
% fs: the sampling frequency [Hz]
% labels: a label for every segment of 30 seconds, for training and testing
% unit: the unit in which the EEG signal is measured
load('./Data/sleepStaging.mat')
load('./Data/EEGcoordsystem.mat');


% Convert the structure array into a cell array
eegCellArray = struct2cell(eeg);
labelsCellArray = struct2cell(labels);

% Use array indexing to access the first element of the cell array (i.e. the "train" element)
traineeg = eegCellArray{1};
trainlabels = labelsCellArray{1};

% Use array indexing to access the second element of the cell array (i.e. the "test" element)
testeeg = eegCellArray{2};
testlabels = labelsCellArray{2};


%% 2.1 Segmentation
% Cut the EEG of both the training and test data in 30s segments.


% Segment the training data into smaller fragments
trainSegments = segmentSignal(traineeg);

% Segment the test data into smaller fragments
testSegments = segmentSignal(testeeg);


% Count the number of occurrences of each label in the trainlabels array
trainCounts = histcounts(trainlabels, {'Wake', 'S2', 'S3'});

% Extract the count for each label from the counts array
trainWakeCount = trainCounts(1);
trainS2Count = trainCounts(2);
trainS3Count = trainCounts(3);

% Calculate the total number of labels in the trainlabels array
trainTotalCount = numel(trainlabels);

% Calculate the percentage of each label in the trainlabels array
trainWakePercent = 100 * trainWakeCount / trainTotalCount;
trainS2Percent = 100 * trainS2Count / trainTotalCount;
trainS3Percent = 100 * trainS3Count / trainTotalCount;

% Count the number of occurrences of each label in the testlabels array
testCounts = histcounts(testlabels, {'Wake', 'S2', 'S3'});

% Extract the count for each label from the counts array
testWakeCount = testCounts(1);
testS2Count = testCounts(2);
testS3Count = testCounts(3);

% Calculate the total number of labels in the testlabels array
testTotalCount = numel(testlabels);

% Calculate the percentage of each label in the testlabels array
testWakePercent = 100 * testWakeCount / testTotalCount;
testS2Percent = 100 * testS2Count / testTotalCount;
testS3Percent = 100 * testS3Count / testTotalCount;

% Create a table with the train percentages
trainTable = array2table([trainWakePercent, trainS2Percent, trainS3Percent], 'VariableNames', {'Wake', 'S2', 'S3'});

% Create a table with the test percentages
testTable = array2table([testWakePercent, testS2Percent, testS3Percent], 'VariableNames', {'Wake', 'S2', 'S3'});

% Display the tables
disp('Train percentages:');
disp(trainTable);
disp('Test percentages:');
disp(testTable);

%% 2.2 Data exploration
% Plot an EEG segment when awake, in sleeping stage 2 (with K-complexes) and
% in sleeping stage 3

% Select a wake segment
wakeSegment = trainSegments{1}

% Select an S2 segment
S2Segment = trainSegments{24}

% Select an S3 segment
S3Segment = trainSegments{44}

figure('Name','Data exploration')
subplot(1,3,1)
% Wake segment
eegPlot(wakeSegment, channels);
title('Wake segment');

subplot(1,3,2)
% S2 segment with at least 1 K-complex
eegPlot(S2Segment, channels);
title('S2 segment');

subplot(1,3,3)
% S3 segment
eegPlot(S3Segment, channels);
title('S3 segment');

%% 2.3 K-complex extraction
% Split the EEG in independent sources using fastica.
% !!!Use the EEG BEFORE segmentation to compute the ICA sources!!!


% Set the seed for the random number generator
rng(1);

[ICs, A, W] = fastica(traineeg);
[ICs_test, A_test, W_test] = fastica(testeeg);


% Segment the estimated sources in segments of 30s each using the same
% function constructed in Section 2.1 

ICSegments = segmentSignal(ICs);
ICSegments_test = segmentSignal(ICs_test);


% Plot the estimated sources of a segment containing K-complexes using
% eegPlot
figure('Name','Estimated sources - temporal')
% Plot of the estimated sources in the time domain

% Select an S2 segment
ICS2Segment = ICSegments{24}


subplot(1,3,2)

% S2 segment with at least 1 K-complex
eegPlot(ICS2Segment, {'IC 1', 'IC 2', 'IC 3', 'IC 4', 'IC 5', 'IC 6', 'IC 7'});
title('IC 24th segment (S2)');

% Create a topoplot of the estimated sources using topoPlot
figure('Name','Estimated sources - spatial')
% Topoplot of the estimated sources

topoPlot(A)

% Select the ICA source that contains the K-complexes most prominently
% hint: the squeeze() function removes irrelevant dimensions of size 1 from
% a given high-dimensional array. After selection, you should have an array
% with two dimensions which only contains 30s segments of the source with 
% K-complexes






%% K-Complex Detection Using Template Matching (training data)

% Load the template
load('./Data/K-complex.mat','template');

% Define the IC source
ic_source = 7;

% Define the segment ranges
segment_ranges = 1:121;  % include all the 121 segments

% Initialize the vector for the number of K-complexes detected per segment with template filtering
K_complexes_per_segment_template = zeros(1, 121);

% Loop through all the segments in the defined ranges
for seg_idx = segment_ranges
    % Select the segment
    ICS2Segment = ICSegments{seg_idx};

    % Select the specified IC source (i.e., the row of the matrix corresponding to the IC source)
    sig = ICS2Segment(ic_source,:)';

    % Normalize the template
    template_norm = template / max(abs(template));

    % Normalize the signal
    sig_norm = sig / max(abs(sig));

    % Compute the cross-correlation
    ccf = xcorr(template, sig, 7500);
    ccf_norm = ccf / max(abs(ccf));
    ccf_norm_cut = ccf_norm(1:7500);

    % Set the threshold value for the cross-correlation
    threshold = 0.999999;

    % Find the sample indices where the cross-correlation crosses the threshold
    onsets = find(ccf_norm_cut >= threshold);

    % Initialize an empty array to store the onsets of the K-complexes
    K_Complex_onsets = [];

    % Iterate over the onsets array
    for i = 1:length(onsets)
        % Check if the current element is the first element or if the difference between the current element and the previous element is greater than 1
        if i == 1 || onsets(i) - onsets(i-1) > 1
            % If either condition is true, add the current element to the K_Complex_onsets array
            K_Complex_onsets = [K_Complex_onsets, onsets(i)];
        end
    end

    % Store the number of K-complexes detected in this segment in the K_complexes_per_segment vector
    K_complexes_per_segment_template(seg_idx) = length(K_Complex_onsets);
end


%% K-Complex Detection Using Matched Filter (training data)

% Load the template
load('./Data/K-complex.mat','template');

% Define the IC source
ic_source = 7;

% Define the segment ranges
segment_ranges = 1:121;  % include all the 121 segments

% Initialize the vector for the number of K-complexes detected per segment with matched filtering
K_complexes_per_segment_matched = zeros(1, 121);

% Loop through all the segments in the defined ranges
for seg_idx = segment_ranges
    % Select the segment
    ICS2Segment = ICSegments{seg_idx};

    % Select the specified IC source (i.e., the row of the matrix corresponding to the IC source)
    sig = ICS2Segment(ic_source,:)';

    % Normalize the template
    template_norm = template / max(abs(template));

    % Normalize the signal
    sig_norm = sig / max(abs(sig));

    % Apply the matched filter
    convolved_signal = fftfilt(template, sig);

    % Set the threshold value for the convolved signal
    threshold = 0.99999999;

    % Find the sample indices where the convolved signal crosses the threshold
    convolved_signal_norm = convolved_signal / max(abs(convolved_signal));
    onsets = find(convolved_signal_norm >= threshold);

    % Initialize an empty array to store the onsets of the K-complexes
    K_Complex_onsets = [];

    % Iterate over the onsets array
    for i = 1:length(onsets)
        % Check if the current element is the first element or if the difference between the current element and the previous element is greater than 1
        if i == 1 || onsets(i) - onsets(i-1) > 1
            % If either condition is true, add the current element to the K_Complex_onsets array
            K_Complex_onsets = [K_Complex_onsets, onsets(i)];
        end
    end

    % Store the number of K-complexes detected in this segment in the K_complexes_per_segment_matched vector
    K_complexes_per_segment_matched(seg_idx) = length(K_Complex_onsets);
end


%% K-Complex Detection Using Template Matching (test data)

% Define the IC source
ic_source = 7;

% Define the segment ranges
segment_ranges = 1:length(ICSegments_test);  % include all the segments

% Initialize the vector for the number of K-complexes detected per segment with template filtering
K_complexes_per_segment_template_test = zeros(1, length(ICSegments_test));

% Loop through all the segments in the defined ranges
for seg_idx = segment_ranges
    % Select the segment
    testSegment = ICSegments_test{seg_idx};

    % Select the specified IC source (i.e., the row of the matrix corresponding to the IC source)
    sig = testSegment(ic_source,:)';

    % Normalize the template
    template_norm = template / max(abs(template));

    % Normalize the signal
    sig_norm = sig / max(abs(sig));

    % Compute the cross-correlation
    ccf = xcorr(template, sig, 7500);
    ccf_norm = ccf / max(abs(ccf));
    ccf_norm_cut = ccf_norm(1:7500);

    % Set the threshold value for the cross-correlation
    threshold = 0.999999;

    % Find the sample indices where the cross-correlation crosses the threshold
    onsets = find(ccf_norm_cut >= threshold);

    % Initialize an empty array to store the onsets of the K-complexes
    K_Complex_onsets = [];

    % Iterate over the onsets array
    for i = 1:length(onsets)
        % Check if the current element is the first element or if the difference between the current element and the previous element is greater than 1
        if i == 1 || onsets(i) - onsets(i-1) > 1
            % If either condition is true, add the current element to the K_Complex_onsets array
            K_Complex_onsets = [K_Complex_onsets, onsets(i)];
        end
    end

    % Store the number of K-complexes detected in this segment in the K_complexes_per_segment vector
    K_complexes_per_segment_template_test(seg_idx) = length(K_Complex_onsets);
end


%% K-Complex Detection Using Matched Filter (test data)

% Define the IC source
ic_source = 7;

% Define the segment ranges
segment_ranges = 1:length(ICSegments_test);  % include all the segments

% Initialize the vector for the number of K-complexes detected per segment with matched filtering
K_complexes_per_segment_matched_test = zeros(1, length(ICSegments_test));

% Loop through all the segments in the defined ranges
for seg_idx = segment_ranges
    % Select the segment
    testSegment = ICSegments_test{seg_idx};

    % Select the specified IC source (i.e., the row of the matrix corresponding to the IC source)
    sig = testSegment(ic_source,:)';

    % Normalize the template
    template_norm = template / max(abs(template));

    % Normalize the signal
    sig_norm = sig / max(abs(sig));

    % Apply the matched filter
    convolved_signal = fftfilt(template, sig);

    % Set the threshold value for the convolved signal
    threshold = 0.99999999;

    % Find the sample indices where the convolved signal crosses the threshold
    convolved_signal_norm = convolved_signal / max(abs(convolved_signal));
    onsets = find(convolved_signal_norm >= threshold);

    % Initialize an empty array to store the onsets of the K-complexes
    K_Complex_onsets = [];

    % Iterate over the onsets array
    for i = 1:length(onsets)
        % Check if the current element is the first element or if the difference between the current element and the previous element is greater than 1
        if i == 1 || onsets(i) - onsets(i-1) > 1
            % If either condition is true, add the current element to the K_Complex_onsets array
            K_Complex_onsets = [K_Complex_onsets, onsets(i)];
        end
    end

    % Store the number of K-complexes detected in this segment in the K_complexes_per_segment_matched_test vector
    K_complexes_per_segment_matched_test(seg_idx) = length(K_Complex_onsets);
end



%% Visualization of K-Complex Detection with Template and Matched filter

% Load the template
load('./Data/K-complex.mat','template');

% Select the 24th segment
ICS2Segment = ICSegments{24};

% Select the third component (i.e., the third row of the matrix)
sig = ICS2Segment(7,:)';
original_sig = ICS2Segment(7,:);  % Copy of the original ICA source
ccf = xcorr(flip(template), sig, 7500);

% Normalize the signals
sig_norm = sig / max(abs(sig));
original_sig_norm = original_sig / max(abs(original_sig));  % Normalize original ICA source
ccf_norm = ccf / max(abs(ccf));
ccf_norm_cut = ccf_norm(1:7500);

% Set the threshold value for the cross-correlated signal
threshold = 0.97;

% Find the sample indices where the cross-correlated signal crosses the threshold
onsets = find(flip(ccf_norm_cut) >= threshold);

% Initialize an empty array to store the onsets of the K-complexes
K_Complex_onsets = [];

% Iterate over the onsets array
for i = 1:length(onsets)
    % Check if the current element is the first element or if the difference between the current element and the previous element is greater than 1
    if i == 1 || onsets(i) - onsets(i-1) > 1
        % If either condition is true, add the current element to the K_complex_onsets array
        K_Complex_onsets = [K_Complex_onsets, onsets(i)];
    end
end

% Recalculate the convolved signal for the selected segment and IC
convolved_signal = conv(sig, template, 'same');

% --- From the peaks in the matched filter's output, compute the onsets of
% the K-complex via a thresholding procedure.
% Set the threshold value for the convolved signal
threshold = 0.999;

% Find the sample indices where the convolved signal crosses the threshold
convolved_signal_norm = convolved_signal / max(abs(convolved_signal));
onsets_conv = find(convolved_signal_norm >= threshold);

% Initialize an empty array to store the onsets of the K-complex
K_Complex_onsets_conv = [];

% Iterate over the onsets array
for i = 1:length(onsets_conv)
  % Check if the current element is the first element or if the difference between the current element and the previous element is greater than 1
  if i == 1 || onsets_conv(i) - onsets_conv(i-1) > 1
    % If either condition is true, add the current element to the K_complex_onsets array
    K_Complex_onsets_conv = [K_Complex_onsets_conv, onsets_conv(i)];
  end
end

% Define sampling frequency
fs = 250; % Hz

% Create a time vector
t = (0:length(sig_norm)-1)/fs; 

% Plotting 
figure('Name','K-complex detection')

% ICA Source plot
subplot(3,1,1)
plot(t, sig_norm,'LineWidth',1.2)
xlim([9 12]) % set x-axis limits from 0 to maximum time
xlabel('Time (s)') % label the x-axis with time in seconds
title('ICA Source')

% Template Filter plot
subplot(3,1,2)
plot(t, sig_norm,'LineWidth',1.2)
hold on
plot(t, flip(ccf_norm_cut))
% Add the detected K-complex onsets to the plot
stem_height = 0.5;
K_Complex_onsets_time = K_Complex_onsets/fs; % convert the K-complex onset indices to time
stem(K_Complex_onsets_time, stem_height * ones(size(K_Complex_onsets)), 'filled', 'k');
legend('Signal', 'Cross-correlation', 'K-complex onsets');
xlim([9 12]) % set x-axis limits from 0 to maximum time
xlabel('Time (s)') % label the x-axis with time in seconds
title('Template Filter')

% Matched Filter plot
subplot(3,1,3)
plot(t, sig_norm);
hold on;
plot(t, convolved_signal_norm);
% Add the detected K-complex onsets to the plot
K_Complex_onsets_time_conv = K_Complex_onsets_conv/fs; % convert the K-complex onset indices to time
stem(K_Complex_onsets_time_conv, stem_height * ones(size(K_Complex_onsets_conv)), 'filled', 'k');
legend('Signal', 'Matched filter output', 'K-complex onsets');
xlim([9 12]) % set x-axis limits from 0 to maximum time
xlabel('Time (s)') % label the x-axis with time in seconds
title('Matched Filter');



%% 2.5 Band power
% Select the f3 channel.

trainF3 = traineeg(4,:);
trainF3Segments = segmentChannel(trainF3);

% Compute the average band power for each segment and each of the 10 bands.

% Define the edges of the ten bands
bands = [0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40];

% Get the number of segments
numSegments = numel(trainF3Segments);

% Preallocate a matrix to store the average power in dB for each segment and each band
avgPower_dB = zeros(numSegments, numel(bands)-1);

% Loop through the segments
for i = 1:numSegments
  % Get the current segment
  trainF3segment = trainF3Segments{i};

  % Compute the PSD of the segment using pwelch
  [PSD, freq] = pwelch(trainF3segment, [], [], 7499.5*2, 250); % specify 7500 frequency bins

  % Convert the PSD values from linear scale to dB scale
  PSD_dB = 10*log10(PSD + 1e-9);
  

  % Loop through the bands and compute the average power in each band
  for j = 1:numel(bands)-1
    % Select the indices of the PSD values corresponding to the current band
    idx = freq >= bands(j) & freq < bands(j+1);

    % Compute the average power in the current band
    avgPower_dB(i,j) = mean(PSD_dB(idx));
    
  end
end

% Transpose the matrix to make the segments the columns of the matrix
avgPower_dB = transpose(avgPower_dB);
%% 2.6 Feature analysis
% Create twelve box plots, showing the statistics of each feature for each
% sleep stage.
figure('Name','Training features')

% Calculate the statistical measures for each feature within each class
% Select the first average band power



kmf_wakeData = K_complexes_per_segment_matched(trainlabels == "Wake");
kmf_S2Data = K_complexes_per_segment_matched(trainlabels == "S2");
kmf_S3Data = K_complexes_per_segment_matched(trainlabels == "S3");
ktf_wakeData = K_complexes_per_segment_template(trainlabels == "Wake");
ktf_S2Data = K_complexes_per_segment_template(trainlabels == "S2");
ktf_S3Data = K_complexes_per_segment_template(trainlabels == "S3");



% Find the length of the longest group
maxLen1 = max(length(kmf_wakeData), max(length(kmf_S2Data), length(kmf_S3Data)));
maxLen2 = max(length(ktf_wakeData), max(length(ktf_S2Data), length(ktf_S3Data)));
maxLen = max(maxLen1, maxLen2);

% Pad the shorter groups with NaN values
kmf_wakeData = padarray(kmf_wakeData, [0, maxLen - length(kmf_wakeData)], NaN, 'post');
kmf_S2Data = padarray(kmf_S2Data, [0, maxLen - length(kmf_S2Data)], NaN, 'post');
kmf_S3Data = padarray(kmf_S3Data, [0, maxLen - length(kmf_S3Data)], NaN, 'post');
ktf_wakeData = padarray(ktf_wakeData, [0, maxLen - length(ktf_wakeData)], NaN, 'post');
ktf_S2Data = padarray(ktf_S2Data, [0, maxLen - length(ktf_S2Data)], NaN, 'post');
ktf_S3Data = padarray(ktf_S3Data, [0, maxLen - length(ktf_S3Data)], NaN, 'post');

% Create a subplot for kmatchedfiltertest
subplot(6, 2, 1)
% Create a boxplot showing the distribution of kmatchedfiltertest within each class
boxplot([kmf_wakeData', kmf_S2Data', kmf_S3Data'],'Labels', {'Wake', 'S2', 'S3'}, 'Symbol', 'o')
% Set the title to kmatchedfiltertest
title('kcountmatchedfilter (training data)')

% Create a subplot for ktemplatefiltertest
subplot(6, 2, 2)
% Create a boxplot showing the distribution of ktemplatefiltertest within each class
boxplot([ktf_wakeData', ktf_S2Data', ktf_S3Data'],'Labels', {'Wake', 'S2', 'S3'}, 'Symbol', 'o')
% Set the title to ktemplatefiltertest
title('kcounttemplatefilter (training data)')

for i = 3:12 % iterate over the 10 average band powers, starting at the 3rd subplot
    avgBandPower = avgPower_dB(i-2,:); % select the (i-2)th average band power
    wakeData = avgBandPower(trainlabels == "Wake");
    S2Data = avgBandPower(trainlabels == "S2");
    S3Data = avgBandPower(trainlabels == "S3");
    
    % Find the length of the longest group
    maxLen = max(length(wakeData), max(length(S2Data), length(S3Data)));
    
    % Pad the shorter groups with NaN values
    wakeData = padarray(wakeData, [0, maxLen - length(wakeData)], NaN, 'post');
    S2Data = padarray(S2Data, [0, maxLen - length(S2Data)], NaN, 'post');
    S3Data = padarray(S3Data, [0, maxLen - length(S3Data)], NaN, 'post');

    % Create a subplot for the (i-2)th average band power
    subplot(6, 2, i)
    % Create a cell array of labels for each column in the input matrix
    labels = cellstr(num2str((1:size(wakeData,2))'));
    
    % Create a boxplot showing the distribution of the average band power within each class
    boxplot([wakeData', S2Data', S3Data'],'Labels', {'Wake', 'S2', 'S3'}, 'Symbol', 'o')
    % Set the title to the (i-2)th average band power
    title(sprintf('Average Band Power %d (training data)', i-2))
    
end




% 2.7.1 Extract the features for each of the test segments.

testF3 = testeeg(4,:);
testF3Segments = segmentChannel(testF3);

% Compute the average band power for each segment and each of the 10 bands.

% Define the edges of the ten bands
bands = [0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40];

% Get the number of segments
testnumSegments = numel(testF3Segments);

% Preallocate a matrix to store the average power in dB for each segment and each band
testavgPower_dB = zeros(testnumSegments, numel(bands)-1);

% Loop through the segments
for i = 1:testnumSegments
  % Get the current segment
  testF3segment = testF3Segments{i};

  % Compute the PSD of the segment using pwelch
  [testPSD, testfreq] = pwelch(testF3segment, [], [], 7499.5*2, 250); % specify 7500 frequency bins

  % Convert the PSD values from linear scale to dB scale
  testPSD_dB = 10*log10(testPSD + 1e-9);

  % Loop through the bands and compute the average power in each band
  for j = 1:numel(bands)-1
    % Select the indices of the PSD values corresponding to the current band
    testidx = testfreq >= bands(j) & testfreq < bands(j+1);

    % Compute the average power in the current band
    testavgPower_dB(i,j) = mean(testPSD_dB(testidx));
    
  end
end

% Transpose the matrix to make the segments the columns of the matrix
testavgPower_dB = transpose(testavgPower_dB);


% Create twelve box plots, showing the statistics of each feature for each
% sleep stage.
figure('Name','Test features')

testkmatchedfiltertest = K_complexes_per_segment_matched_test; % Vector containing K-complexes per segment for matched filter test
testktemplatefiltertest = K_complexes_per_segment_template_test; % Vector containing K-complexes per segment for template filter test


tkmf_wakeData = testkmatchedfiltertest(trainlabels == "Wake");
tkmf_S2Data = testkmatchedfiltertest(trainlabels == "S2");
tkmf_S3Data = testkmatchedfiltertest(trainlabels == "S3");

tktf_wakeData = testktemplatefiltertest(trainlabels == "Wake");
tktf_S2Data = testktemplatefiltertest(trainlabels == "S2");
tktf_S3Data = testktemplatefiltertest(trainlabels == "S3");

% Find the length of the longest group
tmaxLen1 = max(length(tkmf_wakeData), max(length(tkmf_S2Data), length(tkmf_S3Data)));
tmaxLen2 = max(length(tktf_wakeData), max(length(tktf_S2Data), length(tktf_S3Data)));
tmaxLen3 = max(tmaxLen1, tmaxLen2);

% Pad the shorter groups with NaN values
tkmf_wakeData = padarray(tkmf_wakeData, [0, tmaxLen1 - length(tkmf_wakeData)], NaN, 'post');
tkmf_S2Data = padarray(tkmf_S2Data, [0, tmaxLen1 - length(tkmf_S2Data)], NaN, 'post');
tkmf_S3Data = padarray(tkmf_S3Data, [0, tmaxLen1 - length(tkmf_S3Data)], NaN, 'post');

tktf_wakeData = padarray(tktf_wakeData, [0, tmaxLen2 - length(tktf_wakeData)], NaN, 'post');
tktf_S2Data = padarray(tktf_S2Data, [0, tmaxLen2 - length(tktf_S2Data)], NaN, 'post');
tktf_S3Data = padarray(tktf_S3Data, [0, tmaxLen2 - length(tktf_S3Data)], NaN, 'post');

% Create a subplot for kmatchedfiltertest
subplot(6, 2, 1)
% Create a boxplot showing the distribution of kmatchedfiltertest within each class
boxplot([tkmf_wakeData', tkmf_S2Data', tkmf_S3Data'],'Labels', {'Wake', 'S2', 'S3'}, 'Symbol', 'o')
% Set the title to kmatchedfiltertest
title('kcountmatchedfilter (test data)')

% Create a subplot for ktemplatefiltertest
subplot(6, 2, 2)
% Create a boxplot showing the distribution of ktemplatefiltertest within each class
boxplot([tktf_wakeData', tktf_S2Data', tktf_S3Data'],'Labels', {'Wake', 'S2', 'S3'}, 'Symbol', 'o')
% Set the title to ktemplatefiltertest
title('kcounttemplatefilter (test data)')

for i = 3:12 % iterate over the 10 average band powers, starting at the 3rd subplot
    testavgBandPower = testavgPower_dB(i-2,:); % select the (i-2)th average band power
    twakeData = testavgBandPower(testlabels == "Wake");
    tS2Data = testavgBandPower(testlabels == "S2");
    tS3Data = testavgBandPower(testlabels == "S3");
    
    % Find the length of the longest group
    tmaxLen = max(length(twakeData), max(length(tS2Data), length(tS3Data)));
    
    % Pad the shorter groups with NaN values
    twakeData = padarray(twakeData, [0, tmaxLen - length(twakeData)], NaN, 'post');
    tS2Data = padarray(tS2Data, [0, tmaxLen - length(tS2Data)], NaN, 'post');
    tS3Data = padarray(tS3Data, [0, tmaxLen - length(tS3Data)], NaN, 'post');

    % Create a subplot for the (i-2)th average band power
    subplot(6, 2, i)
    % Create a cell array of labels for each column in the input matrix
    tlabels = cellstr(num2str((1:size(twakeData,2))'));
    
    % Create a boxplot showing the distribution of the average band power within each class
    boxplot([twakeData', tS2Data', tS3Data'],'Labels', {'Wake', 'S2', 'S3'}, 'Symbol', 'o')
    % Set the title to the (i-2)th average band power
    title(sprintf('Average Band Power %d (test data)', i-2))
end

%% 2.7.3 Full-model classification

% 2.7.2 Normalization
% Is it required? If yes, normalize the features

% Concatenate the k-complex count row vectors and the average band power matrix along the rows
trainfeatures = vertcat(K_complexes_per_segment_matched, K_complexes_per_segment_template, avgPower_dB);


% Normalize the feature matrix using z-score normalization
trainfeatures_normalized = zscore(trainfeatures);

% Concatenate the testkmatchedfiltertest and testktemplatefiltertest row vectors and the testavgPower_dB matrix along the columns
testfeatures = vertcat(testkmatchedfiltertest, testktemplatefiltertest, testavgPower_dB);

% Normalize the testfeatures matrix using z-score normalization
testfeatures_normalized = zscore(testfeatures);



% Implement a neural network with 2 hidden layers and 10 nodes per hidden
% layer. Train the model on your training features and test it on the test
% features. 

% Convert the categorical array into a matrix of dummy variables
trainlabels_dummy = transpose(dummyvar(trainlabels));

% Create a neural network with two hidden layers and 10 nodes per hidden layer
hiddenLayerSize = [10 10];
outputLayerSize = 3;
net = fitnet([hiddenLayerSize outputLayerSize]);

% Train the network on the normalized training data
[net, tr] = train(net, trainfeatures_normalized, trainlabels_dummy);

% evaluation of the network

% Evaluate the network on the normalized test data
% Each row represents predicted class probability, values represent the
% probability that segment belongs to corresponding class

predictions = net(testfeatures_normalized);
predictions_index = vec2ind(predictions);
predictions_categorical = categorical(predictions_index,[1 2 3], {'S2' 'S3' 'Wake'});


% Calculate the confusion matrix
[cm, order] = confusionmat(testlabels, predictions_categorical);

% Calculate the classification accuracy
accuracy = sum(diag(cm)) / sum(cm(:));

% Calculate sensitivity (true positive rate)
sensitivity = cm(2,2) / (cm(2,2) + cm(2,1));

% Calculate specificity (true negative rate)
specificity = (cm(1,1) + cm(3,3)) / (cm(1,1) + cm(1,2) + cm(3,3) + cm(3,2));

% Create a table of the results
results = table(sensitivity, accuracy, specificity, 'VariableNames', {'Sensitivity', 'Accuracy', 'Specificity'});

% Display the table
disp(results)

figure
% Plot the confusion matrix for the reduced model
plotconfusion(testlabels', predictions_categorical);

% Set the title for the full model
title('Confusion Matrix full model');

%%
% 2.7.4 Feature selection
% Select up to 4 features using the boxplots

% Implement again a neural network with 2 hidden layers and 10 nodes per 
% hidden layer. Train and test it using only the 4 features you selected.
% Average band power 1 & 2 and 7 & 8 distinguish the sleep fases well

% extract average band power 1,2,7 & 8 for training data
trainfeatures_reduced=avgPower_dB([1 2 7 8],:)

% Normalize the feature matrix using z-score normalization
trainfeatures_reduced_normalized = zscore(trainfeatures_reduced);

% extract average band power 1,2,7 & 8 for test data
testfeatures_reduced=testavgPower_dB([1 2 7 8],:)

% Normalize the testfeatures matrix using z-score normalization
testfeatures_reduced_normalized = zscore(testfeatures_reduced);

% Convert the categorical array into a matrix of dummy variables
trainlabels_dummy = transpose(dummyvar(trainlabels));

% Create a neural network with two hidden layers and 10 nodes per hidden layer
hiddenLayerSize = [10 10];
outputLayerSize = 3;
net_reduced = fitnet([hiddenLayerSize outputLayerSize]);

% Train the network on the normalized training data
[net_reduced, tr] = train(net_reduced, trainfeatures_reduced_normalized, trainlabels_dummy);

predictions_reduced = net_reduced(testfeatures_reduced_normalized);
predictions_reduced_index = vec2ind(predictions_reduced);
predictions_reduced_categorical = categorical(predictions_reduced_index,[1 2 3], {'S2' 'S3' 'Wake'});


% Compute the average accuracy and the confusion matrix of this model.
% Calculate the confusion matrix
[cm_reduced, order_reduced] = confusionmat(testlabels, predictions_reduced_categorical);

% Calculate the classification accuracy
accuracy_reduced = sum(diag(cm_reduced)) / sum(cm_reduced(:));

% Calculate sensitivity (true positive rate)
% sensitivity is drastically reduced in comparison with full model
sensitivity_reduced = cm_reduced(2,2) / (cm_reduced(2,2) + cm_reduced(2,1));

% Calculate specificity (true negative rate) 
specificity_reduced = (cm_reduced(1,1) + cm_reduced(3,3)) / (cm_reduced(1,1) + cm_reduced(1,2) + cm_reduced(3,3) + cm_reduced(3,2));

% Create a table of the results
results_reduced = table(sensitivity_reduced, accuracy_reduced, specificity_reduced, 'VariableNames', {'Sensitivity', 'Accuracy', 'Specificity'});


% Show a table with the average accuracy and plot the confusion
% matrices of each model. 

% Display the table
disp(results_reduced)

figure
% Plot the confusion matrix for the reduced model
plotconfusion(testlabels', predictions_reduced_categorical);

% Set the title for the second plot
title('Confusion Matrix reduced model');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BONUS : use the techniques learnt in this course to improve the feature %
% selection. After the selection, compute the average accuracy and        %
% confusion matrix of this new model.                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% BONUS: LDA technique to reduce computational complexity

% Fit LDA model to the normalized training feature matrix
ldaModel = fitcdiscr(trainfeatures_normalized', trainlabels, 'DiscrimType','linear');

% Extract the linear discriminants from the trained model
ldaCoeffs = ldaModel.Coeffs(1,2).Linear;

% Here, we sort the absolute values of the coefficients in descending order
% and select the indices of the 4 largest to identify the most important features
[~, idx] = sort(abs(ldaCoeffs), 'descend');
selected_features_idx = idx(1:4);

% Use only the selected features for training and testing
trainfeatures_reduced = trainfeatures_normalized(selected_features_idx,:);
testfeatures_reduced = testfeatures_normalized(selected_features_idx,:);

% Normalize the feature matrix using z-score normalization
trainfeatures_reduced_normalized = zscore(trainfeatures_reduced')';
testfeatures_reduced_normalized = zscore(testfeatures_reduced');

% Convert the categorical array into a matrix of dummy variables
trainlabels_dummy = dummyvar(trainlabels)';

% Create a neural network with two hidden layers and 10 nodes per hidden layer
hiddenLayerSize = [10 10];
outputLayerSize = 3;  % Adjusted to match the number of classes
net_reduced = fitnet([hiddenLayerSize outputLayerSize]);

% Train the network on the normalized training data
[net_reduced, tr] = train(net_reduced, trainfeatures_reduced_normalized, trainlabels_dummy);

% Make predictions on the test data
predictions_reduced = net_reduced(testfeatures_reduced_normalized');
predictions_reduced_index = vec2ind(predictions_reduced);
predictions_reduced_categorical = categorical(predictions_reduced_index,[1 2 3], {'S2' 'S3' 'Wake'});

% Compute the average accuracy and the confusion matrix of this model.
% Calculate the confusion matrix
[cm_reduced, order_reduced] = confusionmat(testlabels, predictions_reduced_categorical);

% Calculate the classification accuracy
accuracy_reduced = sum(diag(cm_reduced)) / sum(cm_reduced(:));

% Calculate sensitivity (true positive rate) for class 'S3'
sensitivity_reduced = cm_reduced(2,2) / (cm_reduced(2,2) + cm_reduced(2,1));

% Calculate specificity (true negative rate)
% Here, we calculate specificity as the rate of correct classification of 'non-S3' classes
specificity_reduced = (cm_reduced(1,1) + cm_reduced(3,3)) / (cm_reduced(1,1) + cm_reduced(1,2) + cm_reduced(3,1) + cm_reduced(3,3));

% Create a table of the results
results_reduced = table(sensitivity_reduced, accuracy_reduced, specificity_reduced, 'VariableNames', {'Sensitivity', 'Accuracy', 'Specificity'});

% Display the table
disp(results_reduced)

figure
% Plot the confusion matrix for the reduced model
plotconfusion(testlabels', predictions_reduced_categorical);

% Set the title for the second plot
title('Confusion Matrix reduced LDA model');



%% segmentation functions


function segments = segmentSignal(eeg)
  % Calculate the number of samples in each segment
  segmentLength = 250 * 30;

  % Calculate the number of cells needed to divide the signal into segments
  numCells = ceil(size(eeg, 2) / segmentLength);

  % Segment the signal into a cell array, with each cell containing a 7-channel segment
  segments = mat2cell(eeg, 7, segmentLength * ones(1, numCells));
end

function segments = segmentSignal2(eeg)
  % Calculate the number of samples in each segment
  segmentLength = 250 * 30;

  % Calculate the number of cells needed to divide the signal into segments
  numCells = ceil(size(eeg, 2) / segmentLength);

  % Segment the signal into a cell array, with each cell containing a 7-channel segment
  segments = mat2cell(eeg, 5, segmentLength * ones(1, numCells));
end

function segments = segmentChannel(eeg)
  % Calculate the number of samples in each segment
  segmentLength = 250 * 30;

  % Calculate the number of cells needed to divide the signal into segments
  numCells = ceil(size(eeg, 2) / segmentLength);

  % Pad the input matrix with zeros so that its size is a multiple of the segment length
  eeg = [eeg, zeros(1, numCells*segmentLength - size(eeg,2))];

  % Segment the signal into a cell array, with each cell containing a 1-channel segment
  segments = mat2cell(eeg, 1, segmentLength * ones(1, numCells));
end


