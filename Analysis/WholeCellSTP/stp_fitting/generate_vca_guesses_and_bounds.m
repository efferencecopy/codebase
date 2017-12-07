function [guesses, ub, lb] = generate_vca_guesses_and_bounds(model)

check_vca_model(model); % will throw an error if problems

% setup default param values
default_guess_d = 0.6;  % needs to be between [0, 1]
default_guess_d_tau = 0.800; % units are in seconds for the taus
default_guess_f = 1.5;
default_guess_f_tau = 2;

default_ub_d = 1;  % upper bound constrained to be 1
default_ub_d_tau = 7; % units are in seconds for the taus
default_ub_f = 5; 
default_ub_f_tau = 7; 

default_lb_d = 0;  % upper bound constrained to be 1
default_lb_d_tau = 0.050; % units are in seconds for the taus
default_lb_f = 0; 
default_lb_f_tau = 0.050;

N_d_terms = numel(regexpi(model, 'd'));
N_f_terms = numel(regexpi(model, 'f'));

% set up the default guesses
d_guess = ones(1, N_d_terms) .* default_guess_d; 
d_tau_guess = ones(1, N_d_terms) .* default_guess_d_tau; 
f_guess = ones(1, N_f_terms) .* default_guess_f;
f_tau_guess = ones(1, N_f_terms) .* default_guess_f_tau;
guesses = cat(2, d_guess, d_tau_guess, f_guess, f_tau_guess);

% set up the defualt upper bounds
d_ub = ones(1, N_d_terms) .* default_ub_d;
d_tau_ub = ones(1, N_d_terms) .* default_ub_d_tau; % units are in seconds for the taus
f_ub = ones(1, N_f_terms) .* default_ub_f;
f_tau_ub = ones(1, N_f_terms) .* default_ub_f_tau;
ub = cat(2, d_ub, d_tau_ub, f_ub, f_tau_ub);

% set up the lower bounds.
d_lb = ones(1, N_d_terms) .* default_lb_d;
d_tau_lb = ones(1, N_d_terms) .* default_lb_d_tau; % units are in seconds for the taus
f_lb = ones(1, N_f_terms) .* default_lb_f;
f_tau_lb = ones(1, N_f_terms) .* default_lb_f_tau;
lb = cat(2, d_lb, d_tau_lb, f_lb, f_tau_lb);


