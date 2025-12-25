function [hrv_nl] = hrv_nonlinear_edit(nni)

load Nonlinear_defaults

% Define input
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('nni', @(x) isnumeric(x) && ~isscalar(x));
p.addParameter('mse_max_scale', DEFAULT_MSE_MAX_SCALE, @(x) isnumeric(x) && isscalar(x));
p.addParameter('mse_metrics', DEFAULT_MSE_METRICS, @(x) islogical(x) && isscalar(x));
p.addParameter('sampen_r', DEFAULT_SAMPEN_R, @(x) isnumeric(x) && isscalar(x));
p.addParameter('sampen_m', DEFAULT_SAMPEN_M, @(x) isnumeric(x) && isscalar(x));

% Get input
p.parse(nni);
mse_max_scale = p.Results.mse_max_scale;
mse_metrics = p.Results.mse_metrics;
sampen_r = p.Results.sampen_r;
sampen_m = p.Results.sampen_m;

% Create output table
hrv_nl = table;
hrv_nl.Properties.Description = 'Nonlinear HRV Metrics';

%% Preprocess

% Calculate zero-based interval time axis
nni = nni(:);
tnn = [0; cumsum(nni(1:end-1))];

%% Poincare plot

[sd1, sd2] = poincare(nni, 'plot', false);
hrv_nl.SD1 = sd1 * 1000;
hrv_nl.Properties.VariableUnits{'SD1'} = 'ms';
hrv_nl.Properties.VariableDescriptions{'SD1'} = 'NN interval standard deviation along the perpendicular to the line-of-identity';

hrv_nl.SD2 = sd2 * 1000;
hrv_nl.Properties.VariableUnits{'SD2'} = 'ms';
hrv_nl.Properties.VariableDescriptions{'SD2'} = 'NN interval standard deviation along the line-of-identity';

%% DFA-based Nonlinear metrics (short and long-term scaling exponents, alpha1 & alpha2)

% Calcualte DFA
[~, ~, alpha1, alpha2] = dfa_edit(tnn, nni);

% Save the scaling exponents
hrv_nl.alpha1 = alpha1;
hrv_nl.Properties.VariableUnits{'alpha1'} = 'n.u.';
hrv_nl.Properties.VariableDescriptions{'alpha1'} = 'DFA low-scale slope';

hrv_nl.alpha2 = alpha2;
hrv_nl.Properties.VariableUnits{'alpha2'} = 'n.u.';
hrv_nl.Properties.VariableDescriptions{'alpha2'} = 'DFA high-scale slope';

%% Multiscale sample entropy

% Calculate the MSE graph
[ mse_values, scale_axis] = mse_edit(nni, 'mse_max_scale', mse_max_scale, 'sampen_m', sampen_m, 'sampen_r',sampen_r);

% Save the first MSE value (this is the sample entropy).
if ~isempty(mse_values)
    hrv_nl.SampEn = mse_values(1);
    hrv_nl.Properties.VariableUnits{'SampEn'} = 'n.u.';
    hrv_nl.Properties.VariableDescriptions{'SampEn'} = 'Sample entropy';
end

if mse_metrics
    for ii = 1:length(mse_values)
        curr_metric_name = ['MSE' num2str(scale_axis(ii))];
        hrv_nl{:, curr_metric_name} = mse_values(ii);
        hrv_nl.Properties.VariableUnits{curr_metric_name} = 'n.u.';
        hrv_nl.Properties.VariableDescriptions{curr_metric_name} = sprintf('MSE value at scale %d', scale_axis(ii));
    end
end


end




