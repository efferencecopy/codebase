%% initialize the PR655

retval = PR655init()

%% Make a measurement

spd = nan(101, 15);
for i_spd = 1:15
    [spd(:,i_spd), qual, S] = myPR655measspd([],retval);
end


lambda = S(1) : S(2) : (S(3)-1)*S(2) + S(1);

close all
plot(lambda, nanmean(spd,2));

%%

load spd_472_bp; % curtains open
spd_1 = spd;
%spd_1 = mean(spd_1, 2);
spd_1 = bsxfun(@rdivide, spd_1, max(spd_1,[], 1));


load spd_472_lp_bp;
spd_2 = spd;
%spd_2 = mean(spd_2, 2);
spd_2 = bsxfun(@rdivide, spd_2, max(spd_2,[], 1));

load spd_472_lp_bp_hotmirror_1; % curtains closed
spd_3 = spd;
%spd_3 = mean(spd_3, 2);
spd_3 = bsxfun(@rdivide, spd_3, max(spd_3,[], 1));

lambda = S(1) : S(2) : (S(3)-1)*S(2) + S(1);

figure, hold on
plot(lambda, spd_1, 'b')
plot(lambda, spd_2, 'k')
plot(lambda, spd_3, 'g')