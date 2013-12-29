%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
% Initialize a structure (p) that will hold the data for each trial
% We initialize it here.  It gets filled up by repeated calls to 
% ProcessEventList().  Once the "processnowcode" is detected in the 
% data stream, we know we're at the end of a trial so we go through 
% the calculation and plotting routines.  Then we strip the data that
% we just processed from the from the "p" structure.
function p = InitPStruct(processnowcode, headercode)
    if exist('processnowcode', 'var')
        p.processnowcode = processnowcode;
    end
    
    if exist('headercode', 'var')
        p.headercode = headercode;
    end

    p.times = [];
    p.events = [];
    p.processnowflag = 0;
    p.lastprocessed_t = 0;
    p.spikes = {[] [] []};  % Assuming a maximum of three isolated cells
end
