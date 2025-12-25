% Create a time vector for x-axis
timeVector = (1:length(CAP_Epochs(:, 1))) / fs;

figure;
plot(timeVector, CAP_Epochs(:, 1))
ylim([-14, 14])
yticks(-14:2:14);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('EEG Signal');
hold on; ax = gca;

x = CAP_Epochs(:,1);
t = (0:numel(x)-1)/fs;

% Envelope threshold to highlight CAP-like bursts
win   = round(0.5*fs);  % 0.5 s moving window
envx  = movmean(abs(x), win);
madx  = median(abs(envx - median(envx))) + eps;
thr   = median(envx) + 1.5*madx;

mask  = envx > thr;
CC    = bwconncomp(mask);
L     = cellfun(@numel, CC.PixelIdxList);
keep  = find(L >= 1*fs & L <= 60*fs); % 2–60 s bursts

yl = ylim;
for k = 1:numel(keep)
    idx = CC.PixelIdxList{keep(k)};
    x1 = (idx(1)-1)/fs; 
    x2 = (idx(end)-1)/fs;
    patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], ...
          'r','FaceAlpha',0.35,'EdgeColor','none');
end

legend({'EEG','Burst'}, 'Location','northeast');
set(gcf,'Color','w'); grid on;

% Optional: save the annotated figure
print(gcf,'-dpng','-r300','Figure2_CAP_like.png');

eeg_ep = double(CAP_Epochs(:,1));          % 5-min EEG epoch
ecg_ep = double(CAPECG_Epochs(:,1));       % matching ECG epoch
tEEG   = (0:numel(eeg_ep)-1)/fs;

%% --- CAP-like burst shading (same logic as Fig 2, adjustable threshold) ---
win  = round(0.5*fs);                       % 0.5 s smoothing of amplitude envelope
envx = movmean(abs(eeg_ep), win);
madx = median(abs(envx - median(envx))) + eps;
thr  = median(envx) + 2.0*madx;             % make 1.8 or 2.5 if too many/few shades
mask = envx > thr;

CC   = bwconncomp(mask);
L    = cellfun(@numel, CC.PixelIdxList);
keep = find(L >= 1*fs & L <= 60*fs);        % keep spans 2–60 s

cap_windows = zeros(numel(keep),2);
for k = 1:numel(keep)
    idx = CC.PixelIdxList{keep(k)};
    cap_windows(k,:) = [(idx(1)-1)/fs, (idx(end)-1)/fs];   % [start end] in seconds
end

eeg_ep = double(CAP_Epochs(:,1));
ecg_ep = double(CAPECG_Epochs(:,1));
tEEG   = (0:numel(eeg_ep)-1)/fs;

% ---------- CAP-like shading (same logic as Fig 2) ----------
win  = round(0.5*fs);
envx = movmean(abs(eeg_ep), win);
madx = median(abs(envx - median(envx))) + eps;
thr  = median(envx) + 2.0*madx;        % make 1.8 or 2.5 to adjust density
mask = envx > thr;

CC   = bwconncomp(mask);
L    = cellfun(@numel, CC.PixelIdxList);
keep = find(L >= 1*fs & L <= 60*fs);    % 2–60 s spans

cap_windows = zeros(numel(keep),2);
for k = 1:numel(keep)
    I = CC.PixelIdxList{keep(k)};
    cap_windows(k,:) = [(I(1)-1)/fs, (I(end)-1)/fs];
end

fs_ecg = 200;
ecg_signal = CAPECG_Epochs(:,1);   % ECG corresponding to the same epoch
eeg_signal = CAP_Epochs(:,1);

% --- compute simple HRV timeline (instantaneous HR) ---
[qrs_amp_raw, qrs_i_raw, ~] = pan_tompkin(double(ecg_signal), fs_ecg, 0);
RR = diff(qrs_i_raw) / fs_ecg;             % RR intervals (seconds)
tRR = qrs_i_raw(2:end) / fs_ecg;           % time axis for RR intervals

good = RR > 0.50 & RR < 1.50;
RR = RR(good); tRR = tRR(good);
[RR_clean,~] = hampel(RR, 7);
% Remove spikes (ectopy / noise)
RR_med = movmedian(RR, 9, 'omitnan');
bad = abs(RR - RR_med) > 0.20 .* RR_med;
RR_clean = RR;
RR_clean(bad) = RR_med(bad);
HR = 60 ./ RR_clean;
t1Hz = ceil(tRR(1)):1:floor(tRR(end));
HR_1Hz = interp1(tRR, HR, t1Hz, 'pchip');
% smooth with 7-s moving mean (thin but stable)
HR_1Hz = movmean(HR_1Hz, 7, 'omitnan');

% clamp to a plausible sleep range
HR_1Hz = min(max(HR_1Hz, 55), 105);

% --- For plotting on EEG time axis ---
HRi = nan(size(tEEG));
ix = tEEG >= t1Hz(1) & tEEG <= t1Hz(end);
HRi(ix) = interp1(t1Hz, HR_1Hz, tEEG(ix), 'linear');
% --- define CAP-like mask (same as Fig 2) ---
envx = movmean(abs(eeg_signal), round(0.5*fs));
madx = median(abs(envx - median(envx))) + eps;
thr = median(envx) + 1.5*madx;
mask = envx > thr;
CC = bwconncomp(mask);
L  = cellfun(@numel, CC.PixelIdxList);
keep = find(L >= 1*fs & L <= 120*fs);

% --- plotting synchronisation ---
figure('Color','w');
tEEG = (0:numel(eeg_signal)-1)/fs;

subplot(2,1,1)
plot(tEEG, eeg_signal, 'b');
hold on
yl = ylim;
for k = 1:numel(keep)
    idx = CC.PixelIdxList{keep(k)};
    x1 = (idx(1)-1)/fs; x2 = (idx(end)-1)/fs;
    patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], ...
          'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end
ylabel('EEG Amplitude (µV)');
title('CAP Events (red) and EEG signal');

subplot(2,1,2)
plot(tRR, HR, 'k', 'LineWidth', 1);
hold on
yl = ylim;
for k = 1:numel(keep)
    idx = CC.PixelIdxList{keep(k)};
    x1 = (idx(1)-1)/fs; x2 = (idx(end)-1)/fs;
    patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], ...
          'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end
xlabel('Time (seconds)');
ylabel('Heart Rate (bpm)');
title('HRV changes aligned with CAP bursts');
linkaxes(findall(gcf,'type','axes'),'x');

% high-res save
print(gcf,'-dtiff','-r600','Figure3_CAP_HRV_Timeline.tif');
