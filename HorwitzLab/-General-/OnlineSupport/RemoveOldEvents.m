% Function for getting rid of any data that preceeded the 
% "startcode" drop, which should be the first code that
% the online data analysis program cares about.
% It's a little redundant with CleanUpEvents(), but I 
% couldn't see an elegant way to integrate them.
function p = RemoveOldEvents(p, startcode)

    p.processnowflag = 0;
    pncodeidx = find(p.events == startcode,1,'last');
    if (~isempty(pncodeidx))
        p.lastprocessed_t = p.times(pncodeidx);
        L = p.times < p.lastprocessed_t;
        p.times(L) = [];
        p.events(L) = [];
        for i = 1:length(p.spikes)
            spiketimevect = p.spikes{i};
            p.spikes{i} = [spiketimevect(spiketimevect > p.lastprocessed_t)];
        end
    end
end

