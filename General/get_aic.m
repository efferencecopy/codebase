function aic_c = get_aic(raw, fit, num_params)

SSE = sum((raw-fit).^2);
N = numel(raw);
aic = N .* log(SSE./N) + (2 .* num_params+1);
aic_c = aic + ((2.*num_params.*(num_params+1)) ./ (N-num_params-1));