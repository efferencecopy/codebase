function cod = get_r2(raw, fit, num_params)

if ~exist('num_params', 'var') || isempty(num_params)
    norm_fact = 1;
else
    num_observations = numel(raw);
    norm_fact = (num_observations - 1) ./ (num_observations - num_params);
end

SSE = sum((raw-fit).^2);
SST = sum((raw - mean(raw)).^2);

cod = 1 - (norm_fact .* (SSE./SST));

