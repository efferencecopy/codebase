%% WHICH MICE SHOULD CONTRIBUTE TO THE POPULATION ANALYSIS?

% clear out the workspace
fin

% in = {Mouse Name, Site}

in = {'CH_141215_B', 1;...
      'CH_141215_B', 2;...
      'CH_141215_E', 1;...
      'CH_141215_E', 2;...
      'CH_150105_A', 1;...
      'CH_150105_C', 1;...
      'CH_150105_D', 2};



%% LOOP THOUGH EACH MOUSE AND CREATE THE NECESSARY RAW DATA TRACES


% grab the fiber volley pop excel workbook
fname = [GL_DOCUPATH, 'Other_workbooks', filesep, 'fiberVolleyCellList.xlsx'];
[~,txt, raw] = xlsread(fname);
raw(size(txt,1)+1:end, :) = [];
raw(:,size(txt,2)+1:end) = [];
clear txt

% do the analysis
Nexpts = size(in,1);
dat = {};
for i_ex = 1:Nexpts;
    
    % figure out what rows in the work book to pay attention to
    l_mouse = cellfun(@(x) ~isempty(x), regexp(raw(:,1), in{i_ex,1}));
    l_site = [false ; cell2mat(raw(2:end,2)) == in{i_ex,2}];
    l_expt = l_mouse & l_site;
    
    % run the analysis
    [dat{i_ex}, info{i_ex}] = fiberVolleyAnalysis(l_expt, raw, false);
    
    % enter a few other useful things into the structure
    info{i_ex}.mouse =  in{i_ex,1};
    info{i_ex}.opsin = unique(raw(l_expt,8));
    info{i_ex}.ignoreChans = unique(cell2mat(raw(l_expt,6:7)), 'rows');
end



%% PULL OUT THE STATS FOR FURTHER ANALYSIS

% PEAK TO PEAK AMP
% INTEGRAL OF THE TRACE
% save the snippets for plotting

prePulseTime = 0.001; % in sec
postPulseTime = 0.009; % in sec

for i_ex = 1:Nexpts
    
    TF_fields = fieldnames(dat{i_ex});
    Ntfs = numel(TF_fields);
    for i_tf = 1:Ntfs
        
        sampRate = info{i_ex}.(TF_fields{i_tf}).sampRate;
        prePulseSamps = ceil(prePulseTime .* sampRate);
        postPulseSamps = ceil(postPulseTime .* sampRate);
        Nsamps = prePulseSamps + postPulseSamps + 1;
        
        conds = {'FV_Na', 'nbqx_apv_cd2_ttx'};
        for i_cond = 1:numel(conds)
            
            Npulses = sum(info{i_ex}.(TF_fields{i_tf}).pulseOn_idx);
            pulseOn_idx = find(info{i_ex}.(TF_fields{i_tf}).pulseOn_idx);
            dat{i_ex}.(TF_fields{i_tf}).snips.(conds{i_cond}) = {nan(Npulses, Nsamps), nan(Npulses, Nsamps)};
            for i_pulse = 1:Npulses
                
                snip_idx = pulseOn_idx(i_pulse)-prePulseSamps : 1 : pulseOn_idx(i_pulse)+postPulseSamps;
                
               for i_ch = 1:2;
                    
                    % deal with some exceptions
                    if ~info{i_ex}.ignoreChans(i_ch)
                        continue
                    end
                    
                    % a strange case where HS1 was set to Vclamp, and so
                    % HS2's data got put in the first column...
                    if strcmpi(info{i_ex}.mouse, 'CH_150105_D')
                        if i_ch == 1; error('something went wrong'); end
                        i_ch = 1;
                    end
                            
                    
                    % pull out the snippet and store it for each pulse in
                    % the train
                    snippet_full = dat{i_ex}.(TF_fields{i_tf}).(conds{i_cond})(snip_idx ,i_ch);
                    dat{i_ex}.(TF_fields{i_tf}).snips.(conds{i_cond}){i_ch}(i_pulse,:) = snippet_full;
                    
                    
                    % figure out the peak to peak amplitude of the FV_na or
                    % the diff from baseline for the _ttx case;
                    snippet_pulse = snippet_full(prePulseSamps+1 : end);
                    snippet_baseline = snippet_full(1:prePulseSamps);
                    
                    switch conds{i_cond}
                        case 'FV_Na'
                            [~, troughidx] = min(snippet_pulse);
                            [~, peakidx] = max(snippet_pulse);
                            assert(peakidx > troughidx, 'ERROR: negativity does not lead the positivity')
                            %assert(peakidx./sampRate < 0.006, 'ERROR: positivity occurs too late');
                            
                            trough = mean(snippet_pulse(troughidx-2:troughidx+2));
                            peak = mean(snippet_pulse(peakidx-2:peakidx+2));
                            
                            pk2tr = peak-trough;
                            
                            dat{i_ex}.(TF_fields{i_tf}).peaks.(conds{i_cond}).pk2tr{i_ch}(i_pulse) = pk2tr;
                            
                        case 'nbqx_apv_cd2_ttx'
                            baseline = mean(snippet_baseline);
                            [~, minidx] = min(snippet_pulse);
                            peakval = mean(snippet_pulse(minidx-2:minidx+2));
                            diffval = peakval - baseline;
                            dat{i_ex}.(TF_fields{i_tf}).peaks.(conds{i_cond}).diffval{i_ch}(i_pulse) = diffval;
                            
                    end
                    
                    
                    
                end
            end
            
        end
        
    end
    
end


%% MAKE SOME PLOTS

close all

%cond = 'FV_Na';
cond = 'nbqx_apv_cd2_ttx';

for i_ex = 1:Nexpts
    
    for i_ch = 1:2;
        
        % deal with data and channels that need to be ignored.
        if ~info{i_ex}.ignoreChans(i_ch)
            continue
        else
            % a strange case where HS1 was set to Vclamp, and so
            % HS2's data got put in the first column...
            if strcmpi(info{i_ex}.mouse, 'CH_150105_D')
                if i_ch == 1; error('something went wrong'); end
                i_ch = 1;
            end
        end
        
        figure
        figName = sprintf('%s: site %d, ch %d', info{i_ex}.mouse, in{i_ex,2}, i_ch);
        set(gcf, 'name', figName);
        
        
        
        %
        % plot the raw (Average) traces
        %
        TF_fields = fieldnames(dat{i_ex});
        Ntfs = numel(TF_fields);
        for i_tf = 1:Ntfs
            
            tmp_raw = dat{i_ex}.(TF_fields{i_tf}).snips.(cond){i_ch};
            
            tt = (0:size(tmp_raw,2)-1) ./ info{i_ex}.(TF_fields{i_tf}).sampRate;
            tt = (tt - prePulseTime) ./ 1000;
            
            subplot(2, Ntfs, i_tf)
            cmap = colormap('copper');
            cidx = round(linspace(1, size(cmap,1), size(tmp_raw,1)));
            cmap = cmap(cidx,:);
            set(gca, 'colororder', cmap, 'NextPlot', 'replacechildren');
            
            plot(tt, tmp_raw', 'linewidth', 2)
            axis tight
            title(['TF = ', TF_fields{i_tf}(4:end), 'Hz'])
            xlabel('time (ms)')
            if i_tf==1;
                ylabel('LFP amp')
            end
            
        end
        
        
        
        %
        % plot the summary stat as a function of pulse number
        %
        switch cond
            case 'nbqx_apv_cd2_ttx'
                
                subplot(2,Ntfs, Ntfs+1), hold on,
                for i_tf = 1:Ntfs
                   
                    diffvals = dat{i_ex}.(TF_fields{i_tf}).peaks.(cond).diffval{i_ch};
                    diffvals = abs(diffvals);
                    plot(1:numel(diffvals), diffvals, 'o-', 'color', cmap(i_tf,:), 'linewidth', 2)
                end
                xlabel('Pulse number')
                ylims = get(gca, 'ylim');
                set(gca, 'ylim', [0, ylims(2)]);
                
                
            case 'FV_Na'
                
                subplot(2,Ntfs, Ntfs+1), hold on,
                for i_tf = 1:Ntfs
                   
                    pk2tr = dat{i_ex}.(TF_fields{i_tf}).peaks.(cond).pk2tr{i_ch};
                    plot(1:numel(pk2tr), pk2tr, 'o-', 'color', cmap(i_tf,:), 'linewidth', 2)
                end
                xlabel('Pulse number')
                ylims = get(gca, 'ylim');
                set(gca, 'ylim', [0, ylims(2)]);
        end
        
        
        
    end
end


