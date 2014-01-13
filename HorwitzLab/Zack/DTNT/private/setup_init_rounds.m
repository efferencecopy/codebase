function setup_init_rounds(threshs, color_dirs, trajectories)
global gl
threshs = threshs';
lmsmat = mkbasis(color_dirs')';

% Each row is a color direction for consistency with NeuroThreshOnline.m
% Making entries into the trialspecs array: round 1, (3 color dirs)
for ii = 1:size(lmsmat,1)
    gl.dtnt.ts(ii).colordir = lmsmat(ii,:)./norm(lmsmat(ii,:));
    gl.dtnt.ts(ii).measuredthreshold = threshs(ii);
    gl.dtnt.ts(ii).OOG = threshs(ii) == max(trajectories{ii});
    gl.dtnt.nspecs = gl.dtnt.nspecs + 1;
end

% Making entries into the trialspecs array: round 2, (4 color dirs)
signmat = [(fullfact([2,2])-1.5)*2 ones(4,1)];  % assumes three color dirs
parentOOGs = [gl.dtnt.ts(1:gl.dtnt.nspecs).OOG];
for ii = 1:size(signmat,1)
    gl.dtnt.nspecs = gl.dtnt.nspecs + 1;
    parentvertices = diag(signmat(ii,:).*threshs)*lmsmat;
    v = mean(parentvertices);
    gl.dtnt.ts(gl.dtnt.nspecs).colordir = v./norm(v);
    gl.dtnt.ts(gl.dtnt.nspecs).predictedthreshold = norm(v);
    gl.dtnt.ts(gl.dtnt.nspecs).parentvertices = parentvertices;
    gl.dtnt.ts(gl.dtnt.nspecs).parentOOGs = parentOOGs;
end
