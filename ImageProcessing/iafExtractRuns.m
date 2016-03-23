function out = iafExtractRuns(fnames, relevantTrialAttributes)

% 'fnames' must be an Nx2 cell array. The first column contains the .mat
%          files and the second column contains the .tiff data
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
fprintf('Loading in data for %d runs\n', size(fnames, 1))
fprintf(' Loading the Matlab (trial by trial) data files\n')
Nruns = size(fnames, 1);
for i_run = 1:Nruns;
    
    % Loading/Opening the datasets ('fnames' files)
    info_img{i_run} = loadTIFF_info(fnames{i_run, 2}); % Structure array with one element for each image in the file (per run)
    
    % Loading in the trial attributes. In some cases the .mat file contains
    % multiple structures, so just take the one called 'input'
    tdat = load(fnames{i_run, 1});
    info_run{i_run} = tdat.input; 
    
    % Defining the matrix attributes (# of rows/columns), according to 'info'
    rows(i_run) = info_img{i_run}.height_pix;
    columns(i_run) = info_img{i_run}.width_pix;
    
    % Important parameters...
    Nimage(i_run) = double(info_run{i_run}.counter{end}); % Determining the total # of images (per run)
    info_img{i_run}.Nframes = Nimage(i_run);
    Ntrials(i_run) = numel([info_run{i_run}.counter{:}]); % Determine the number of trials according the MATLAB file
    
    NframesON(i_run) = info_run{i_run}.nScansOn; % Determine the number of frames during stim on
    NframesOFF(i_run) = info_run{i_run}.nScansOff; % Determine the number of frames during the ISI
   
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

% Preallocate space...
totalNtrials = sum(Ntrials);
allON = repmat({nan(rows, columns, NframesON)}, totalNtrials, 1);
allOFF = repmat({nan(rows, columns, NframesOFF)}, totalNtrials, 1);
invalidTrials = false(totalNtrials, 1);


for i_run = 1:Nruns;
    
    % grab the trial frame counter
    all_trial_starts = [0, info_run{i_run}.counter{:}];
    
    % call loadTIFF
    fprintf(' Extracting tiff file %d of %d. ', i_run, Nruns); tic;
    run_dfof = loadTIFF(fnames{i_run,2}, info_img{i_run}); %not actually dfof but setting up for in-place operations
    fprintf('%.1f seconds \n', toc)
    run_dfof = double(run_dfof);
    
    % convert to dfof
    frameRate = 1000 ./ double(info_run{i_run}.frameImagingRateMs);
    window_size_sec = (NframesON + NframesOFF)./frameRate * 2; % two trials worth of seconds
    fprintf('  Now computing dFoF. '); tic
    run_dfof = dfof_from_tiffstack(run_dfof, frameRate, window_size_sec);
    fprintf('%.0f seconds \n', toc)
    
    
    % extract each trial (on and off response separately)
    for i_trial = 1: Ntrials(i_run);
        
        % Define a start and stop position
        idx_start = all_trial_starts(i_trial)+1;
        idx_stop = all_trial_starts(i_trial+1);
        
        % Make sure it's a valid trial
        INVALID = false;
        if idx_stop <= idx_start;
            INVALID = true;
        end
        
        % the the trial data from the run
        frameNums = idx_start : idx_stop;
        Trial_dfof = run_dfof(:, :, frameNums);
        
        
        
        % figure out what the trial number is (cumulative over runs)
        if i_run > 1
            trials_from_previous_runs = cumsum(Ntrials(1 : i_run-1));
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
        
        
        % Define a start and stop position for ON/OFF parts of trial
        frameNumsOFF = 1 : NframesOFF;
        frameNumsON = (NframesOFF + 1) : (NframesOFF + NframesON);
        assert((NframesOFF + NframesON) == size(Trial_dfof, 3), 'ERROR: inconsistent number of frames');
        
        % Extracting OFF parts of trial
        if i_trial == 1
            % do nothing b/c the "off" response for the Nth stimulus comes
            % from the N+1 trial. And the first "off" portion is not in
            % response to any stimulus
        else
            % put the "off" data for the Nth trial in the slot for the N-1
            % trial. This is b/c the "off" response is causally related to
            % what came before it, but the "off" period is in the beginnign
            % of each trial. A consequence is that the last trial has no
            % "off" response.
            allOFF{combined_trial_num-1} = Trial_dfof(:,:,frameNumsOFF);
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