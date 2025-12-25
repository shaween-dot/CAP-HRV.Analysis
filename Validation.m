T = readtable('ML_V.xlsx');
disp(T)

metrics = { ...
    'AVNN', ...
    'SDNN', ...
    'RMSSD', ...
    'LF', ...
    'HF', ...
    'LF_HF', ...
    'TotalPower', ...
    'SampEn', ...
    'DFA1', ...
    'DFA2'};

R2   = zeros(numel(metrics),1);
RMSE = zeros(numel(metrics),1);
MAE  = zeros(numel(metrics),1);
Bias = zeros(numel(metrics),1);

for i = 1:numel(metrics)

    x = T.([metrics{i} '_MATLAB']);      % MATLAB pipeline values
    y = T.([metrics{i} '_PZ']);   % PhysioZoo GUI values
    % Remove any NaNs just in case
    idx = isfinite(x) & isfinite(y);
    x = x(idx);
    y = y(idx);

    % R²
    R  = corr(x, y);
    R2(i) = R^2;

    % RMSE
    RMSE(i) = sqrt(mean((x - y).^2));

    % Optional but recommended
    MAE(i)  = mean(abs(x - y));
    Bias(i) = mean(x - y);
end

ValidationTable = table( ...
    metrics', R2, RMSE, MAE, Bias, ...
    'VariableNames', {'Metric','R2','RMSE','MAE','Bias'});

disp(ValidationTable)



