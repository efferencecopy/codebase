function [colorDirs_out, contrasts_out] = convert_to_Rstar_colordirs(colorDirs_in, contrasts_in, mon)

% the 'mon' input structure must contain the following fields:
%     mon.Mmtx
%     mon.bkgndrgb
%     mon.rgb2Rstar
%     mon.bkgndlms_Rstar


contrasts = cat(1,contrasts_in{:});
nColors = size(colorDirs_in,1);


% initialize the output arguments
Rstar_colorDir = nan(size(colorDirs_in));
Rstar_contrasts = nan(size(contrasts));

% convert all the colors into lms_Rstar color dirs. In the process, figure
% out what the new contrasts should be
for i_clr = 1:nColors
    
    tmp_colorDir_smj = colorDirs_in(i_clr,:);
    tmp_colorDir_smj = tmp_colorDir_smj ./ norm(tmp_colorDir_smj); %make sure it's a unit vector
    
    % scale the colorDir by all the contrasts used. now the colorDirs are
    % in CC units
    tmp_colorDirs_smj = bsxfun(@times, tmp_colorDir_smj,  contrasts(i_clr,:)');
    
    % determine the device specific lms (in smj units)
    bkgndlms_deviceSpecific_smj = mon.Mmtx * mon.bkgndrgb(:);
    gaborlms_deviceSpecific_smj = bsxfun(@times, bkgndlms_deviceSpecific_smj', (1+tmp_colorDirs_smj));
    
    % determine what rgb direction this would correspond to on the
    % moniotor used in the monkey experiments
    gaborrgb_deviceSpecific = mon.Mmtx \ gaborlms_deviceSpecific_smj';
    
    % convert the device specific rgbs into lms units for the cone
    % model
    lms_Rstar = mon.rgb2Rstar * gaborrgb_deviceSpecific;
    diffval = bsxfun(@minus, lms_Rstar, mon.bkgndlms_Rstar);
    tmp_colorDir_Rstar = bsxfun(@rdivide, diffval, mon.bkgndlms_Rstar);
    
    % determine the LMS color dir and the contrasts
    CCs_Rstar = sqrt(sum(tmp_colorDir_Rstar.^2, 1));
    colorDir_Rstar_unitvec = bsxfun(@rdivide, tmp_colorDir_Rstar, CCs_Rstar);
    l_nonzeroContrast = CCs_Rstar > 1e-10;
    maxvals = max(colorDir_Rstar_unitvec(:,l_nonzeroContrast)',[],1);
    minvals = min(colorDir_Rstar_unitvec(:,l_nonzeroContrast)',[],1);
    assert(all([maxvals-minvals] < 1e-10), 'ERROR, something went wrong calculating color dirs')
    
    % store the information
    Rstar_colorDir(i_clr,:) = colorDir_Rstar_unitvec(:,end)';
    Rstar_contrasts(i_clr,:) = CCs_Rstar;
end


% package the outputs back into the 'gab' structure
colorDirs_out = Rstar_colorDir;
contrasts_out = mat2cell(Rstar_contrasts, ones(size(contrasts,1),1), size(contrasts,2));






