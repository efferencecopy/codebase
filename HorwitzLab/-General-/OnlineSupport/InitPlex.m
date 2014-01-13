%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize a Plexon client process
function s = InitPlex()
    s = 0;
    while ~s
        s = PL_InitClient(0);
    end
end