% Contents:
%
% succones: Data from the Schnapf, Nunn, and Baylor paper. Cone action
% spectra measured with suction electrodes
%
% T_cones_synthgh1: A set of cone fundamentals obtained by converting the 
% Schnapf et al. fits to their absorbance data to absorptance (OD = 0.3)
% and adding the Boettner macaque pre-retinal filters. (All media, direct
% measurements).
%
% T_cones_synthgh2: A set of cone fundamentals obtained from the Stockman
% and Sharpe (2000) 10 degree cone fundamentals, with macular pigment
% reduced to 0 and lens multiplied to OD = 1 @ 400 nm . See TenDegFundStuff.m, section 7.2. 
%
% T_cones_synthgh3: A set of cone fundamentals obtained from the Stockman
% and Sharpe (2000) 10 degree cone fundamentals, with macular pigment
% reduced to 0, lens reduced to 0, and the Boettner macaque pre-retinal
% filters applied with OD = 1 @ 400 nm (and shifted to have zero density
% at long wavelengths, consistent with other lens absorption spectra).
% See TenDegFundStuff.m, section 7.2. 
%
% T_cones_synthgh4: A set of cone fundamentals obtained from the Stockman
% and Sharpe (2000) 10 degree cone fundamentals, with macular pigment
% reduced to 0, lens multiplied to OD = 1 @ 400 nm, and S-cone pigment
% spectrum moved to slightly longer wavelength (peak at 430 instead of 
% 425 nm). See TenDegFundStuff.m, section 7.2. (Seems to make zero
% difference)
%
% T_cones_synth5: Like #2 but lens density at 400 nm is only 0.75.
%
% dens_lens_ss: Based on den_lens_ssf, which comes packaged inside the 
% PsychophysicsToolbox. Downsampled to be on the [380:5:780] lattice and 
% absorption at 380 and 385 nm estimated by linear extrapolation on the
% log10 values. See TenDegFundStuff.m, section 7.2.
%
% dens_mac_ss: Based on macss_5.csv, which came from the CVRL.org website.
% Downsampled to be on the [380:5:780] lattice and absorption at 380 and 
% 385nm estimated by linear extrapolation on the log10 values.
% See TenDegFundStuff.m, section 7.2. 
%
% wrattennd1: Transmittance spectrum for a Kodak 1 log-unit neutral density
% filter. ZLB acquired these data on 8/13/12 from page 182 of Kodak filter
% book. Three variables in this file.
%      wratten_nd1_dig: digitized from graph (wls, %trans)
%      wratten_nd1_trans: transcribed from text (wls, %trans)
%      wratten_nd1: interpolating filter spectrum to be on the usual 
%               380:5:780 lattice (see TenDegFundStuff.m, section 11)