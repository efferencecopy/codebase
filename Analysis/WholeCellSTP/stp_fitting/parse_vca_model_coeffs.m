function [d, tau_d, f, tau_f] = parse_vca_model_coeffs(params_vec, model)

check_vca_model(model); % will bonk if there are problems

N_d_terms = numel(regexpi(model, 'd'));
N_f_terms = numel(regexpi(model, 'f'));

start_idx = 1;

stop_idx = N_d_terms;
d = params_vec(start_idx : stop_idx);
start_idx = start_idx + N_d_terms;

stop_idx = start_idx + N_d_terms - 1;
tau_d = params_vec(start_idx : stop_idx);
start_idx = start_idx + N_d_terms;

stop_idx = start_idx + N_f_terms - 1;
f = params_vec(start_idx : stop_idx);
start_idx = start_idx + N_f_terms;

stop_idx = start_idx + N_f_terms - 1;
tau_f = params_vec(start_idx : stop_idx);

assert(stop_idx == numel(params_vec), 'ERROR: may have iterated incorrectly through the params vector')