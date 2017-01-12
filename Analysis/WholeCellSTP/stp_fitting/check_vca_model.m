function is_correct = check_vca_model(model)

assert(ischar(model), 'ERROR: model must be specified as a string')
assert(numel(model)>0, 'ERROR: model must specify at least one ''d'' or ''f'' term.')

% check to see if the model only contains Ds and Fs
is_correct = isempty(regexpi(model, '[^df]'));
assert(is_correct, 'ERROR: model string may only contain the letters ''d'' and ''f''.')