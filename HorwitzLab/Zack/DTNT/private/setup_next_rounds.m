function setup_next_rounds(threshs, color_dirs, trajectories, filename)
global gl
threshs = threshs';
lmsmat = mkbasis(color_dirs')';

LINPREDTOL = 0.3;

prev_color_dirs = reshape([gl.dtnt.ts(1:gl.dtnt.nspecs).colordir], 3, gl.dtnt.nspecs)';
for ii = 1:size(lmsmat,1)
    L = logical(softEq((lmsmat(ii,:)./norm(lmsmat(ii,:)))*prev_color_dirs', 1, 10));

    if ~any(L)
        warning('DTNTOnline:CantFindColorDir', ...
            sprintf('Cannot find color direction in the trial specs: (% .4f,% .4f,% .4f)', lmsmat(ii,:))); %#ok<SPWRN>
        continue
    elseif sum(L) > 1
        error('DTNTOnline:MultipleColorDirs', 'Multiple identical color directions detected.');
    elseif ~isempty(gl.dtnt.ts(L).measuredthreshold)
        error('DTNTOnline:MultipleThresholds', 'Already measured a threshold in this direction');
    end
    gl.dtnt.ts(L).measuredthreshold = threshs(ii);
    gl.dtnt.ts(L).OOG = threshs(ii) == max(trajectories{ii});

    threshratio = gl.dtnt.ts(L).measuredthreshold ./ gl.dtnt.ts(L).predictedthreshold;
    if abs(log(threshratio)) > abs(log(1+LINPREDTOL)) % making new nodes (directions to search)
        grandparentvertices = gl.dtnt.ts(L).parentvertices;
        grandparentOOGs = gl.dtnt.ts(L).parentOOGs;

        for jj = 1:3
            parentOOGs = grandparentOOGs;
            parentOOGs(jj) = gl.dtnt.ts(L).OOG;
            if ~all(parentOOGs)
                gl.dtnt.nspecs = gl.dtnt.nspecs + 1;
                parentvertices = grandparentvertices;
                parentvertices(jj,:) = gl.dtnt.ts(L).colordir .* gl.dtnt.ts(L).measuredthreshold;
                v = mean(parentvertices);
                gl.dtnt.ts(gl.dtnt.nspecs).colordir = v./norm(v);
                gl.dtnt.ts(gl.dtnt.nspecs).predictedthreshold = norm(v);
                gl.dtnt.ts(gl.dtnt.nspecs).parentvertices = parentvertices;
                gl.dtnt.ts(gl.dtnt.nspecs).parentOOGs = parentOOGs;
                gl.dtnt.ts(gl.dtnt.nspecs).measuredthreshold = [];
            end
        end
    else
        fprintf('%s: Threshold in direction (% .4f,% .4f,% .4f) is consistent with a linear prediction\n', ...
            filename, lmsmat(ii,:));
    end
    if abs(log(threshratio)) > log(3)
        warning('DTNTOnline:AberrantThresh', ...
            sprintf('Aberrant threshold point: %s along (% .4f,% .4f,% .4f)', filename, gl.dtnt.ts(L).colordir)); %#ok<SPWRN>
    end
end
