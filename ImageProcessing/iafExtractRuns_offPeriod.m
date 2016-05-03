function out = iafExtractRuns_offPeriod(name_mat, name_img, relevantTrialAttributes, DECORR)

% 
%
% 'relevantTrialAttributes' must be a cell array of strings denoting the 
%          trial attributes that should be used to distinguish trial types.

% initalize the outputs
out.uniqueTrlTypes = [];
out.trlTypeHeader = [];
out.OnResponse = {};
out.OffResponse = {};

% initalize some local variables
trialAttributes = [];

% load the .mat files to get trial by trial data
fprintf('Loading in data for %d runs\n', numel(name_mat))
fprintf(' Loading the Matlab (trial by trial) data files\n')
Nruns = numel(name_mat);
for i_run = 1:Nruns;
    
    % Loading/Opening the datasets ('fnames' files)
    info_img{i_run} = loadTIFF_info(name_img{i_run}); % Structure array with one element for each image in the file (per run)
    
    % Loading in the trial attributes. In some cases the .mat file contains
    % multiple structures, so just take the one called 'input'
    tdat = load(name_mat{i_run});
    info_run{i_run} = tdat.input; 
    
    % Defining the matrix attributes (# of rows/columns), according to 'info'
    rows(i_run) = info_img{i_run}.height_pix;
    columns(i_run) = info_img{i_run}.width_pix;
    
    % Important parameters...
    Nimage(i_run) = double(info_run{i_run}.counter{end}); % Determining the total # of images (per run)
    info_img{i_run}.Nframes = Nimage(i_run);
    Ntrials(i_run) = numel([info_run{i_run}.counter{:}]); % Determine the number of trials according the MATLAB file
    
    NframesON(i_run) = double(info_run{i_run}.nScansOn); % Determine the number of frames during stim on
    NframesOFF(i_run) = double(info_run{i_run}.nScansOff); % Determine the number of frames during the ISI
   
end

% make sure things are consistent across runs.
assert(numel(unique(rows))==1, 'ERROR: rows are not equal across runs')
assert(numel(unique(columns))==1, 'ERROR: columns are not equal across runs')
assert(numel(unique(NframesON))==1, 'ERROR: NframesON are not equal across runs')
assert(numel(unique(NframesOFF))==1, 'ERROR: NframesOFF are not equal across runs')
rows = rows(1);
columns = columns(1);
NframesON = NframesON(1);
NframesOFF = NframesOFF(1);
frameNumsON = 1 : NframesON;
frameNumsOFF = (NframesON + 1) : (NframesON + NframesOFF);

% other useful info for now and later analyses
frameRate = 1000 ./ double(info_run{i_run}.frameImagingRateMs);
out.frameRate = frameRate;

% Preallocate space...
totalNtrials = sum(Ntrials);
allON = repmat({nan(rows, columns, NframesON)}, totalNtrials, 1);
allOFF = repmat({nan(rows, columns, NframesOFF)}, totalNtrials, 1);
invalidTrials = false(totalNtrials, 1);

for i_run = 1:Nruns;
    
    % call loadTIFF
    fprintf(' Extracting tiff file %d of %d. ', i_run, Nruns); tic;
    img_raw = loadTIFF(name_img{i_run}, info_img{i_run});
    fprintf('%.1f seconds \n', toc)
    img_raw = double(img_raw);
    
    % grab the trial frame counter
    all_trial_starts = double([0, info_run{i_run}.counter{:}]);
    
    % extract each trial (on and off response separately)
    for i_trial = 1: Ntrials(i_run);
        
        % Define a start and stop position
        idx_start = all_trial_starts(i_trial)+1+NframesOFF;
        if i_trial < Ntrials(i_run)
            idx_stop = all_trial_starts(i_trial+1)+NframesOFF;
        else
            idx_stop = all_trial_starts(i_trial+1);
        end
            
        % Make sure it's a valid trial
        INVALID = false;
        if idx_start >= idx_stop;
            INVALID = true;
            error('Need to deal with this case')
        end
        
        % grab the trial data from the run and convert to dfof units
        frameNums = idx_start : idx_stop;
        Trial_dfof = dfof_from_trial(img_raw, frameNums, frameRate);
        if DECORR
            Trial_dfof = decor_mean_subtract(Trial_dfof);
        end
        
        
        % figure out what the trial number is (cumulative over runs)
        if i_run > 1
            trials_from_previous_runs = sum(Ntrials(1 : i_run-1));
        else
            trials_from_previous_runs = 0;
        end
        trials_from_current_run = i_trial;
        combined_trial_num = trials_from_previous_runs + trials_from_current_run;
        
        
        % grab the trial attributes
        for i_attrib = 1:numel(relevantTrialAttributes)
            fieldName = relevantTrialAttributes{i_attrib};
            trialAttributes(combined_trial_num, i_attrib) = double(info_run{i_run}.(fieldName){i_trial});
        end
        
        % so that you can cull invalid trials after importing all the data
        if INVALID
            invalidTrials(combined_trial_num) = true;
        end
        
        % error checking:
        if i_trial < Ntrials
            assert((NframesOFF + NframesON) == size(Trial_dfof, 3), 'ERROR: inconsistent number of frames');
        else
            assert(NframesON == size(Trial_dfof, 3), 'ERROR: inconsistent number of frames');
        end
        
        
        
        % Extracting OFF parts of trial
        if i_trial < Ntrials(i_run)
            % The last trial does not have a true off response, so just
            % aggregate the off responses for trials 1 to Ntrials-1
            allOFF{combined_trial_num} = Trial_dfof(:,:,frameNumsOFF);
        end
        
        % Extracting ON parts of trial
        allON{combined_trial_num} = Trial_dfof(:,:,frameNumsON); 
    end
end

assert(~any(invalidTrials), 'ERROR: invalid trials were found, but not deleted');

% Defining trial types
out.uniqueTrlTypes = unique(trialAttributes, 'rows');
out.trlTypeHeader = relevantTrialAttributes;
N_types = size(out.uniqueTrlTypes, 1);

%Separating Trials by trial types (t_types)
for i_t_types = 1 : N_types;
    
    tType_list = ismember(trialAttributes, out.uniqueTrlTypes(i_t_types,:), 'rows');
    assert(size(unique(trialAttributes(tType_list,:),'rows'),1)==1, 'ERROR: mismatched trial types?')
    
    out.OffResponse{i_t_types} = allOFF(tType_list);
    out.OnResponse{i_t_types} = allON(tType_list);
end

fprintf('Done importing data\n');

end


function dfof = dfof_from_trial(img_raw, trialFrameNums, frameRate)
    
    % use the last few seconds of the preceeding OFF period to convert F
    % into dfof
    Nframes_bkgnd = ceil(frameRate .* 5); % 5 seconds worth of bkgnd time
    bkgnd_idx = [(trialFrameNums(1)-Nframes_bkgnd) : (trialFrameNums(1)-1)];
    zeroFrames = sum(sum(img_raw(:,:,bkgnd_idx),1),2) == 0;
    assert(~any(zeroFrames), 'ERROR: found some background (Fo) frames that were zero')
    img_bkgnd = mean(img_raw(:,:,bkgnd_idx), 3);
    
    dfof = bsxfun(@minus, img_raw(:,:,trialFrameNums(1):trialFrameNums(end)), img_bkgnd);
    dfof = bsxfun(@rdivide, dfof, img_bkgnd);
end

function dfof = decor_mean_subtract(dfof)
    % crazy idea: subtract the mean across all pixels to eliminate image wide
    % noise that has nothing to do with IAF or hemodynamics. Make sure that the
    % thing you subtract off has the same sigma as each pixel time series.
    xbar = mean(mean(dfof,1),2);
    sigma = std(dfof,[],3);
    scaleFactor = sigma ./ std(xbar(:)); % a scale factor to equate sigma on a pix by pix basis
    rsub = bsxfun(@times, xbar, scaleFactor);
    assert(all(all((sigma - std(rsub,[],3))<1e-10)), 'ERROR: sigmas are not the same')
    
    dfof = dfof - rsub;
    
end




