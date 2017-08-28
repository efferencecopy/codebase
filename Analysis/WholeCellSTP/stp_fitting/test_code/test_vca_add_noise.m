function noisy_fake_data_with_ns = test_vca_add_noise(fake_data_with_ns, fake_data_no_ns, train_noise_prcnt)

sigma = fake_data_no_ns .* train_noise_prcnt;
noise = normrnd(zeros(size(sigma)), sigma);

% don't do anything to the first pulse, since noise has already been added
%noise(:,1) = 0;
warning('Hack to make strictly additive noise')

% add the noise to the fake data
noisy_fake_data_with_ns = fake_data_with_ns + noise;

% make sure the data are strictly  positive
noisy_fake_data_with_ns = max(noisy_fake_data_with_ns, 0);