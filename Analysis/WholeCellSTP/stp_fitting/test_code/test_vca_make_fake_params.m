function fake_params = test_vca_make_fake_params(model)

[~, ub, lb] = generate_vca_guesses_and_bounds(model);
fake_params = unifrnd(lb, ub);

[d, tau_d, f, tau_f] = parse_vca_model_coeffs(fake_params, model);

[d, didx] = sort(d);
tau_d = tau_d(didx);

[f, fidx] = sort(f);
tau_f = tau_f(fidx);

fake_params = [d, tau_d, f, tau_f];