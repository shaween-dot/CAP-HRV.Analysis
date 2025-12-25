% According to https://www.medlink.com/articles/cyclic-alternating-pattern,
% CAP is characterized by phase A: sleep phasic activity and phase B:
% return to EEG background.
% According to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8638906/ phase
% A and B can last between 2 and 60 seconds.
%% Start
clear;clc;close all;
SubjectDetails = readtable('MostafaSDP2 Data-2.xlsx');                          % Read the table with subject details.
ExpectedEvents = SubjectDetails.TotalApneas;
%% Initialize Variables
Subject = 2;
Filenames = string(lower(SubjectDetails.Patient_)); % Study number, used for file names, extracted from table
SubjectEDF = 'pat1.edf';                                      % Setting subject edf file name.
[EDFInfo,EDFData] = edfread(SubjectEDF);                          % Ri =reading and extracting info from the EDF file.
EDFLabels = string(EDFInfo.label);
Fs = EDFInfo.frequency;
fs = Fs(EDFLabels == "EMGMiddle");
FsECG = Fs(contains(EDFLabels,"ECG"));
FsECG = FsECG(1);
FsEEG = fs; %Fs(EDFLabels == "EEG");
SignalLength = length(EDFData(strcmpi(EDFInfo.label,"EMGMiddle"),:));
ECGData = EDFData(strcmpi(EDFInfo.label,"ECG1"),:);
if isempty(ECGData)
    ECGData = EDFData(strcmpi(EDFInfo.label,"ECG1"),:);
end
ECGData2 = EDFData(strcmpi(EDFInfo.label,"ECG"),:);
% Remove tailing zeros in ECG data. Not necessary for EEG data
ECGData2 = ECGData2(1:SignalLength/(fs/FsECG));
ECGData2 = ECGData2 - median(ECGData2);
%%

C3A2 = EDFData(strcmpi(EDFInfo.label, "C3A2"),:);
C4A1 = EDFData(strcmpi(EDFInfo.label, "C4A1"),:);
F3A2 = EDFData(strcmpi(EDFInfo.label, "F3A2"),:);
F4A1 = EDFData(strcmpi(EDFInfo.label, "F4A1"),:);
O1A2 = EDFData(strcmpi(EDFInfo.label, "O1A2"),:);
O2A1 = EDFData(strcmpi(EDFInfo.label, "O2A1"),:);
%% Filter Design
FiltOrder = 4; % Keep filter order minimum for faster live filtering
Notch = fdesign.notch('N,F0,Q', FiltOrder, 50, 50, fs); % Notch Filter @ 50 Hz due to power line interference
Filter1 = design(Notch, 'butter');
BandPass = fdesign.bandpass('N,F3dB1,F3dB2', FiltOrder, 0.5, 30, fs); % Bandpass Filter @ 0.5-30 Hz delta to beta
Filter2 = design(BandPass, 'butter');
[pECG,qECG] = rat(fs/FsECG);
ECGData2 = resample(ECGData2,pECG,qECG,100); % Resample ECG data to same sampling frequency as EEG

% Apply Filter2 to EEG data
C3A2 = filter(Filter2,C3A2);
C4A1 = filter(Filter2,C4A1);
F3A2 = filter(Filter2,F3A2);
F4A1 = filter(Filter2,F4A1);
O1A2 = filter(Filter2,O1A2);
O2A1 = filter(Filter2,O2A1);

% Define the delta and theta bandpass filter designs
DeltaBand = fdesign.bandpass('N,F3dB1,F3dB2', FiltOrder, 0.5, 4, fs); % Delta band (0.5-4 Hz)
ThetaBand = fdesign.bandpass('N,F3dB1,F3dB2', FiltOrder, 4, 8, fs);  % Theta band (4-8 Hz)

% Design the delta and theta bandpass filters
DeltaFilter = design(DeltaBand, 'butter');
ThetaFilter = design(ThetaBand, 'butter');

% Apply the delta and theta filters to EEG data
C3A2_delta = filter(DeltaFilter, C3A2);
C4A1_delta = filter(DeltaFilter, C4A1);
F3A2_delta = filter(DeltaFilter, F3A2);
F4A1_delta = filter(DeltaFilter, F4A1);
O1A2_delta = filter(DeltaFilter, O1A2);
O2A1_delta = filter(DeltaFilter, O2A1);


C4A1_delta = zscore(C4A1_delta);
SignalLengthNew = SignalLength;
epoch_duration = 300*fs;
while mod(SignalLengthNew,epoch_duration) ~= 0
    SignalLengthNew = SignalLengthNew-1;
end
C4A1_small = C4A1_delta(1:SignalLengthNew);
C4A1reshaped = reshape(C4A1_small,epoch_duration,[]);
ECG_small = ECGData2(1:SignalLengthNew);
ECGreshaped = reshape(ECG_small,epoch_duration,[]);
numberOfBursts = zeros(size(C4A1reshaped,2),1);

thresholdFactor = 1/3;  % Adjust this factor as needed
for bb = 1:size(C4A1reshaped, 2)
    % Calculate the threshold based on the median amplitude
    medianAmplitude = median(C4A1reshaped(:, bb));
    threshold = thresholdFactor * medianAmplitude;

    % Identify spans that exceed the threshold
    aboveThreshold = (C4A1reshaped(:, bb) > threshold);

    % Find spans that are long enough
    isLongEnough = bwareafilt(aboveThreshold, [1 * fs, 300 * fs]);

    % Count the number of spans (bursts) that are long enough
    [~, numberOfBursts(bb)] = bwlabel(isLongEnough);
end

% numel(find(numberOfBursts))
CAP_Epochs = C4A1reshaped(:,find(numberOfBursts)); %#ok<FNDSB>
CAPECG_Epochs = ECGreshaped(:,find(numberOfBursts)); %#ok<FNDSB>
% plot(CAP_Epochs(:,2))

% Count the number of CAP events
numCAPEvents = numel(find(numberOfBursts));

disp(['Number of CAP Events: ' num2str(numCAPEvents)]);

% Create a time vector for x-axis
timeVector = (1:length(CAP_Epochs(:, 1))) / fs;

figure;
plot(timeVector, CAP_Epochs(:, 1))
ylim([-14, 14])
yticks(-14:2:14);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('EEG Signal during delta burst');


% Identify non-CAP epochs
nonCAP_indices = setdiff(1:size(ECGreshaped, 2), find(numberOfBursts));

% Create a table for non-CAP epochs
nonCAPTable = table();

for i = nonCAP_indices
    % Extract the ECG data for non-CAP epoch
    nonCAP_ECG_epoch = ECGreshaped(:, i);

    % Create a unique variable name for each non-CAP epoch
    varName = ['NonCAP_ECG_epoch_' num2str(i)];

    % Append the current non-CAP epoch to the nonCAPTable
    nonCAPTable.(varName) = nonCAP_ECG_epoch;
end

% Convert the table to a double array
nonCAP = table2array(nonCAPTable);

% Assuming CAPECG_Epochs is a cell array of signals
all_signals = reshape(CAPECG_Epochs, [], 1);
all_signals2 = reshape(nonCAP, [], 1);

% Display the size of the resulting array
disp(['Size of nonCAPArray: ' num2str(size(nonCAP))]);
totallength = numel(C4A1_delta);
caprate = numel(all_signals)/totallength *100;
disp(['cap rate: ' num2str(caprate)]);

% --- Auto-save entire workspace at the end ---
save('MyWorkspace.mat');     % saves all variables in the current workspace
fprintf('Workspace saved to %s\\MyWorkspace.mat\n', pwd);
