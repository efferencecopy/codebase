%% Building a keepratio -> # stims mapping
load('IsoSamp\isosurf_params.mat');
kr_stims_map = zeros(100, 4, 'uint16');
for sfidx = 3:4
    count = 0;
    plotlims = sedna.plotlims(:,sfidx);
    model = sedna.model(:,sfidx);
    
    n_meshpts = 70;
    [x,y,z] = meshgrid( ...
        linspace(-plotlims(1), plotlims(1), n_meshpts), ...
        linspace(-plotlims(2), plotlims(2), n_meshpts), ...
        linspace(-plotlims(3), plotlims(3), n_meshpts));
    
    xyz = [x(:) y(:) z(:)];
    v = sum(abs(xyz * reshape(model(2:end), [3 3])) .^ model(1), 2);
    fv = isosurface(x, y, z, reshape(v, size(x)), 1);
    
    p0 = min(fv.vertices(:,1:3));
    p1 = max(fv.vertices(:,1:3));
    
    for keepratio = linspace(.1, .0001, 100)
        try
            [~,~,face] = surf2mesh(fv.vertices, fv.faces, p0, p1, ...
                keepratio, 5, [], []);
        catch E
            if strcmp(E.identifier, 'MATLAB:unassignedOutputs') || ...
                    strcmp(E.message, 'node file is missing!')
                logger('No stimuli generated (problems with keep ratio?)!\n');
                break
            else
                rethrow(E);
            end
        end
        count = count + 1;
        kr_stims_map(count,sfidx) = uint16(size(face, 1));
    end
end