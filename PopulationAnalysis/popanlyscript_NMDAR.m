%% IMPORT THE NECESSARY DATA INTO A datDB STRUCTURE

fin

mdb = initMouseDB('new'); % clear out the junk from previous datDB runs.

% buid a structure of params from physiology notes. The globals are used
% during execution of the physiology_notes script
global GL_ADD_TO_MDB GL_SUPPRESS_ANALYSIS
GL_ADD_TO_MDB = true; %#ok<NASGU>
GL_SUPPRESS_ANALYSIS = true; %#ok<NASGU>
physiology_notes % run to store prams in MDB.

% cd to the directory of workbooks and import the spreadsheet of cells
workbookpath = [GL_DOCUPATH, filesep, 'Other_workbooks', filesep, 'Cell_Library'];
workbooksheet = 'NMDAR';
[~, ~, raw] = xlsread(workbookpath, workbooksheet);
header = raw(1,:);
mouseNames = raw(2:end,1);
siteNumber = raw(2:end,2);
goodNeurons(:,1) = cat(1,raw{2:end, 3});
goodNeurons(:,2) = cat(1,raw{2:end, 4});
goodNeurons = logical(goodNeurons);
neuronType = raw(2:end, 6:7);
HVA = raw(2:end, 5);
layer = [raw(2:end, 8); raw(2:end, 9)];

%
% build a population data structure by performing the appropriate analysis.
% I'll use the params.name field to build each new entry into the DB.
% Finish by saving the structure (possibly with a date string)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initalize things that I care about
[dat.NMDAR.raw_pA,...
 dat.NMDAR.raw_mV,...
 dat.ivcurve.mV,...
 dat.ivcurve.mv_corrected,...
 dat.ivcurve.pA,...
 dat.tvec] = deal(repmat({[] []}, numel(mouseNames), 1));

% iterate over the mice in the cell library (some mice get analyzed
% multiple times if there were multiple experiments per mouse)
for ex = 1:numel(mouseNames)
    
    ex_mouseName = mouseNames{ex};
    ex_siteNum = siteNumber{ex};
    mdb = initMouseDB(false, true);
    [~, idx] = mdb.search(ex_mouseName);
    
    params = mdb.mice{idx}.popAnly{ex_siteNum};
    params.fxns = {@anlyMod_optoIV, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
    close all; drawnow
    
    %
    % add the data to an array for the AMPA/NMDA
    % Excitation/Inhibition ratio stuff. Adding 'NMDAR' to the 'group'
    % array allows the script to pull out NMDAR IV curve data
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if isfield(params.ivdat, 'NMDAR')
       
        dat.NMDAR.raw_pA(ex,:) = params.ivdat.NMDAR.raw;
        dat.NMDAR.raw_mV(ex,:) = params.ivdat.NMDAR.vhold;
        dat.tvec{ex,1} = params.ivdat.tvec;
        dat.ivcurve.mV(ex, :) = params.ivdat.NMDAR.ivcurve.mV;
        dat.ivcurve.mv_corrected(ex, :) = params.ivdat.NMDAR.ivcurve.mV_corrected;
        dat.ivcurve.pA(ex, :) = params.ivdat.NMDAR.ivcurve.pA;
    end
    
end % iterate over recording sites

% package all the useful things into a single structure
dat.hva = HVA;
dat.goodNeurons = goodNeurons;
dat.neuronType = neuronType;
dat.mice = mouseNames;
dat.siteNum = cat(1, siteNumber{:});
dat.layer = layer;


% create grouping lists for HVAs
hvas = repmat(dat.hva, 2,1);
hvaList.('pm') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'pm'));
hvaList.('lm') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'lm'));
hvaList.('al') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'al'));
hvaList.('und') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'und'));

%
% create grouping lists for neuron type
neuronType = dat.neuronType(:);
neuronType = cellfun(@num2str, neuronType, 'uniformoutput', false);
l_IN = cellfun(@(x) ~isempty(x), regexpi(neuronType, 'in'));
l_SOM = cellfun(@(x) ~isempty(x), regexpi(neuronType, 'som'));
typeList.IN = l_IN | l_SOM;
typeList.SOM = l_SOM;
typeList.PY = cellfun(@(x) ~isempty(x), regexpi(neuronType, 'py'));
typeList.und = ~typeList.IN & ~typeList.PY;


%
% create a grouping list for layer
layer = dat.layer(:);
layer = cellfun(@num2str, layer, 'uniformoutput', false);
layerList.L_23 = cellfun(@(x) ~isempty(x), regexpi(layer, '2/3'));
layerList.L_4 = cellfun(@(x) ~isempty(x), regexpi(layer, '4'));
layerList.L_5 = cellfun(@(x) ~isempty(x), regexpi(layer, '5'));


% save the data and the grouping lists
originalDir = pwd;
cd(GL_POPDATPATH);
save('popAnly_NMDAR.mat', 'dat', 'hvaList', 'typeList', 'layerList')
cd(originalDir);

% be nice and return these variables to their default values
GL_ADD_TO_MDB = false;
GL_SUPPRESS_ANALYSIS = false;


%% PLOT THE IV CURVES FOR EACH AREA

fin


% load in the pre-saved population data
load([GL_POPDATPATH, 'popAnly_NMDAR.mat'])

%
% PARAMETERS FOR ANALYSIS
%
MINNUMVHOLD = 6; % minimum number of vholds a data set needs to be included
l_23 = layerList.L_23;
l_valid = dat.goodNeurons(:) & typeList.PY;

% create a container for the scaled IV curves
scaled.al.pA = {};
scaled.al.mV = {};
scaled.pm.pA = {};
scaled.pm.mV = {};

% plot the IV curves for PY cells
tmp_mV = dat.ivcurve.mV(:);
tmp_mV_corrected = dat.ivcurve.mv_corrected(:);
tmp_pA = dat.ivcurve.pA(:);
tmp_mice = [dat.mice(:); dat.mice(:)];
tmp_siteNum = [dat.siteNum(:); dat.siteNum(:)];
l_AL = hvaList.al;
l_PM = hvaList.pm;
l_und = hvaList.und;
l_toplot = l_valid;
h_raw = figure; hold on,
h_corrected = figure; hold on,
for a = find(l_toplot)'
   
   % figure out the plot colors
   if l_AL(a)
       [clr_fit, clr_raw] = hvaPlotColor('AL');
   elseif l_PM(a)
       [clr_fit, clr_raw] = hvaPlotColor('PM');
   elseif l_und(a)
       [clr_fit, clr_raw] = hvaPlotColor('und');
   else
       continue % do nothing if this isn't AL or PM
   end
   
   % normalize the IV curve by the value at +50 mV. 
   plt_mV = tmp_mV{a};
   plt_mV_corrected = tmp_mV_corrected{a};
   plt_pA = tmp_pA{a};
   
   if numel(plt_pA) < MINNUMVHOLD
       continue
   end
   
   % plot the "corrected" version first b/c it doesn't need any
   % normalization, and thus all the cells will get plotted (not just the
   % ones that have a Vhold = 50mV.
   figure(h_corrected)
   plot(plt_mV_corrected, plt_pA, 'o-', 'color', clr_raw)
   
   
   % now normalize the raw values.
   idx = softEq(plt_mV, 50, 0);
   if sum(idx) == 1
       normfact_pA = plt_pA(idx);
       
   elseif ~any(idx)
       
       warning('No Vhold = 50 mV found')
       
       % try to interpolate b/w two points...
       template = sign(plt_mV - 50);
       template = conv(template, [-1 1]);
       idx = find(template == -2);
       
       if isempty(idx);
           continue
       else
           xx = plt_mV([idx-1, idx]);
           yy = plt_pA([idx-1, idx]);
           xq = 50;
           normfact_pA = interp1(xx, yy, xq);
           warning('Estimating normalization factor')
       end

   else
       error('unknown number of Vhold=50 conditions')
   end
   plt_pA = plt_pA ./ normfact_pA;
   

   figure(h_raw)
   p = plot(plt_mV, plt_pA, 'o-', 'color', clr_raw);
   printTitle = @(a,b,c) title(sprintf('%s, cell %d',c{1},c{2}));
   set(p, 'buttonDownFcn', {printTitle, {tmp_mice{a}, tmp_siteNum(a)}})
   t = get(get(p, 'parent'), 'title');
   
      
   % store the scaled values for later
   if l_AL(a)
      scaled.al.mV{end+1} = plt_mV;
      scaled.al.pA{end+1} = plt_pA;
   elseif l_PM(a)
      scaled.pm.mV{end+1} = plt_mV;
      scaled.pm.pA{end+1} = plt_pA;
   end
   
   
end


figure(h_raw)
xlabel('Voltage (mV)')
ylabel('Current (pA)')
title('IV curve, raw')
set(t, 'interpreter', 'none');

figure(h_corrected)
xlabel('Voltage')
ylabel('Current (pA)')
title('IV curve, corrected for steady state Rs')



%
% Average scaled IV curve
%

tested_mV = [-80 -60 -40 -20  0 15 17 50];
scaled.al.avg = repmat({[]}, 1, numel(tested_mV));
scaled.pm.avg =  repmat({[]}, 1, numel(tested_mV));
group = {'al', 'pm'};
for i_group = 1:2
    for i_cell = 1:numel(scaled.(group{i_group}).mV)
        mV = scaled.(group{i_group}).mV{i_cell};
        pA = scaled.(group{i_group}).pA{i_cell};
        for i_vhold = 1:numel(mV)
            idx_avg = softEq(mV(i_vhold), tested_mV, 0);
            if sum(idx_avg) == 0
                continue
            elseif sum(idx_avg)>1
                error('too many matches')
            end
            scaled.(group{i_group}).avg{find(idx_avg)}(end+1) = pA(i_vhold);
            
        end
    end
end

figure, hold on,
for i_group = 1:2
     [clr_fit, clr_raw] = hvaPlotColor(group{i_group});
     plot(tested_mV, cellfun(@mean, scaled.(group{i_group}).avg), 'o-', 'color', clr_raw)
end
    


%% NMDA RAW CURRENTS



fin


% load in the pre-saved population data
load([GL_POPDATPATH, 'popAnly_NMDAR.mat'])
l_23 = layerList.L_23;
l_valid = dat.goodNeurons(:) & l_23 & typeList.PY;



tmp_Raw = dat.NMDAR.raw_pA(:);
tmp_tvec = [dat.tvec; dat.tvec];
tmp_names = repmat(dat.mice,2,1);
tmp_mV = dat.NMDAR.raw_mV(:);
Nfigs = ceil(sum(l_valid)/4);
cellList = find(l_valid);
idx = 1;
for ii_fig = 1:Nfigs
    
    figure;
    set(gcf, 'position', [130           5        1177         801])
    
    for ii_plt = 1:4
        if numel(cellList) < idx; continue; end
        if isempty(tmp_Raw{cellList(idx)}); idx = idx+1; continue; end
        
        tvec = tmp_tvec{cellList(idx),1};
        l_window = (tvec>0.0035) & (tvec<0.100);
        
        tmp_pA = tmp_Raw{cellList(idx)};
        tmp_pA = cat(1, tmp_pA{:})';
        tmp_pA(~l_window,:) = [];
        
        % a complicated way of picking the correct normalization value
        % regardless of the polarity of the current
        maxVal = max(tmp_pA,[],1);
        minVal = min(tmp_pA, [],1);
        l_out = abs(maxVal) > abs(minVal);
        l_in = abs(maxVal) < abs(minVal);
        peakVal = nan(size(maxVal));
        peakVal(l_out) = maxVal(l_out);
        peakVal(l_in) = minVal(l_in);
        
        % now do the normalization
        tmp_pA = bsxfun(@rdivide, tmp_pA, peakVal);
        
        
        
        subplot(2,2,ii_plt);
        plot(tvec(l_window), tmp_pA)
        t = title(tmp_names{cellList(idx)});
        set(t, 'interpreter', 'none')
        
        leg = tmp_mV{cellList(idx)};
        %leg = mat2cell(leg, ones(1,size(leg,1)), ones(1,size(leg,2)));
        leg = cellfun(@(x) sprintf('%.1f', x), leg, 'uniformoutput', 0);
        legend(leg)
    
        idx = idx+1;
    end
    
    
end


