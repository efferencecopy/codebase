function fnames = getDetectionFiles(minNumTrials);

fnames = {}; %in case we return early
PARADIGMID = 210;
% First making a gigantic dirstruct
dirstruct = dir;
for i = 1:size(dirstruct,1)
    if (dirstruct(i).isdir && ~strncmp(dirstruct(i).name, '.',1))
        cd(dirstruct(i).name);
        dirstruct = [dirstruct; dir];
        cd('../')
    end
end

counter = 1;
for i = 1:size(dirstruct,1)
    if (length(dirstruct(i).name) < 4)
        continue;
    end

    dirstruct(i).name
    if (dirstruct(i).name(end-3:end) == '.nex') %has to be a nex file
        stro = nex2stro(findfile(dirstruct(i).name, pwd));
        if isfield(stro,'sum') %can't be corrupt
            disp(stro.sum.paradigmID);
            if (stro.sum.paradigmID == PARADIGMID) %must be detection file
                disp(stro.sum.exptParams.expt_meth)
                if stro.sum.exptParams.expt_meth == 1; %must be mocs
                    if any(strcmp('sig001a', stro.sum.rasterCells))
                        cntrstLevInd = strmatch('cntrst_lev', [stro.sum.trialFields(1,:)]);
                        colorDirInd = strmatch('color_dir', [stro.sum.trialFields(1,:)]);
                        colorDirs = unique(stro.trial(:,colorDirInd));
                        numContrasts = length(unique(stro.trial(:,cntrstLevInd)));
                        trialsPerCond = nan(numContrasts, size(colorDirs,1));
                        for clr = 1:size(colorDirs,1)
                            l_color = stro.trial(:,colorDirInd) == clr;
                            for cntrst = 1:numContrasts
                                l_cntrst = stro.trial(:,cntrstLevInd) == cntrst;
                                if cntrst == 1;
                                    trialsPerCond(cntrst, clr) = sum(l_cntrst);
                                else
                                    trialsPerCond(cntrst, clr) = sum(l_cntrst&l_color);
                                end
                            end
                        end

                        %add the file to the list if it has enough trials
                        disp(trialsPerCond)
                        disp(min(min(trialsPerCond)));
                        try
                        if ~isempty(trialsPerCond >= minNumTrials) && all(trialsPerCond(:) >= minNumTrials);
                            fnames{counter} = dirstruct(i).name;
                            counter = counter+1;
                        end
                        catch
                            keyboard
                        end
                    end
                end
            end
        end
    end
end