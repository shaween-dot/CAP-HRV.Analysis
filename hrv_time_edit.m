function [hrv_td] = hrv_time_edit(nni,pnn_thresh_ms)

% Convert to milliseconds
nni = nni * 1000;

hrv_td = table;
hrv_td.Properties.Description = 'Time domain HRV metrics';

hrv_td.AVNN = mean(nni);
hrv_td.Properties.VariableUnits{'AVNN'} = 'ms';
hrv_td.Properties.VariableDescriptions{'AVNN'} = 'Average NN interval duration';

hrv_td.SDNN = sqrt(var(nni));
hrv_td.Properties.VariableUnits{'SDNN'} = 'ms';
hrv_td.Properties.VariableDescriptions{'SDNN'} = 'Standard deviation of NN interval duration';

hrv_td.RMSSD = sqrt(mean(diff(nni).^2));
hrv_td.Properties.VariableUnits{'RMSSD'} = 'ms';
hrv_td.Properties.VariableDescriptions{'RMSSD'} = 'The square root of the mean of the sum of the squares of differences between adjacent NN intervals';

hrv_td.pNNx = 100 * sum(abs(diff(nni)) > pnn_thresh_ms) / (length(nni)-1);
hrv_td.Properties.VariableUnits{'pNNx'} = '%';
hrv_td.Properties.VariableDescriptions{'pNNx'} = sprintf('Percent of NN interval differences greater than %.1fmilliseconds', pnn_thresh_ms);
hrv_td.Properties.VariableNames{'pNNx'} = ['pNN' num2str(floor(pnn_thresh_ms))]; % change name to e.g. pNN50

hrv_td.SEM = hrv_td.SDNN / sqrt(length(nni));
hrv_td.Properties.VariableUnits{'SEM'} = 'ms';
hrv_td.Properties.VariableDescriptions{'SEM'} = 'Standard error of the mean NN interval';
