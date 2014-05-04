function mdb = addPopAnlyParamsToMDB(params)

global GL_DOCUPATH

% remove the axon binary file stuff, because this doesn't need to get saved
% in the mdb.
params = rmfield(params, 'ax');

% open the mdb, figure out which mouse the popParams belong to
mdb = initMouseDB(false, true);
[~, idx] = mdb.search(params.mouse);

% do some simple error checking
assert(strcmpi(mdb.mice{idx}.name, params.mouse), 'ERROR saving to mdb, there was a mismatch in mouse names');

% save the popParams information to the mdb. 
mdb.mice{idx}.popAnly{params.cellNum} = params;
save([GL_DOCUPATH, 'mouseDB.mat'], 'mdb')