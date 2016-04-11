function dat = wcstp_compile_data(exinfo, hidx)
%
% INPUTS:
%
% exinfo:  a cell array of experiment attributes, one attribute per column.
%          Contains things like file names etc.
% 
% hidx:    a structure, that contains the lookup info for the fields in
%          'exinfo'. 
%
% OUTPUTS:
%
% dat.qc.Rs: series resistance for HS1,2 [1xNtrials]
% dat.qc.p1amp: amplitude for first pulse of each sweep for HS1,2 [1xNtrials]
% dat.qc.vhold: measured holding potential for each sweep, HS1,2 [1xNtrials]
% 
% dat.dcsteps.Vm_raw: raw voltage traces during stimulus on period for HS1,2 [Nsweeps x Ntime]
% dat.dcsteps.Icmd: magnitude of current injection on a sweep by sweep basis. [Nsweeps x 1]
% 
% dat.expt.[stimType].raw.snips: pA trace surrounding the EPSC for HS1,2 [Nsweeps x Ntime]
% dat.expt.[stimType].stats.amp: EPSC amplitude for HS1,2 [Nsweeps, 1]
% dat.expt.[stimType].realTrialNum: the actual trial number for each of the sweeps
% dat.expt.[stimType].tdict: the trialtype dictionary entry for this condition
% 
% dat.info.opsin
% dat.info.pretime.vclamp
% dat.info.pretime.dcsteps
% dat.info.posttime.vclamp
% dat.info.posttime.dcsteps
% dat.info.sampRate.vclamp
% dat.info.sampRate.iclamp


keyboard

