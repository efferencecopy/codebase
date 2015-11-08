function stro = blk2stro(nev, varargin)

% syntax for input arguments:
% stro = blk2stro(pathToNEV, 'NSx', pathToNSx', 'trialdef', {trial def modules})

% trial def could be something like: ns4::trl_start

inargs = parseInputs(nev, varargin);
keyboard





% import the nev
nev = openNEV(nev, 'read', 'nosave', 'nomat');

% import the nsX data


% determine the trial start/stop times using the event times (spike data
% from NEV). Do some error checking, and perhaps double check with the
% continuous data in the ns6 data (which is sampled at high rate)


% iterate over the data and build the 'ras' field.  The dimensions of 'ras'
% are <Ntrials x Nchannels>, where each channel can have a spike time
% chanel, a spike WF,  and one or more continuous channels. Deal out each
% data packet to a cell array.
%
% **** need to convert things into uV ******


% build the 'sum' field. There could be a main exptParams field that comes 


% build the 'trial' field. Include the trials start/stop times. Determine
% the pulse: amplitude, width, frequency, ipi, stim on time (relative to
% trl start).


end % main function

function p = parseInputs(nev, optionalInputs)
    p = inputParser;
    p.addRequired('nev',      @(x) (ischar(x) && any(regexpi(x, '.nev'))))
    p.addParameter('ns1', '', @(x) (ischar(x) && any(regexpi(x, '.ns1'))))
    p.addParameter('ns2', '', @(x) (ischar(x) && any(regexpi(x, '.ns2'))))
    p.addParameter('ns3', '', @(x) (ischar(x) && any(regexpi(x, '.ns3'))))
    p.addParameter('ns4', '', @(x) (ischar(x) && any(regexpi(x, '.ns4'))))
    p.addParameter('ns5', '', @(x) (ischar(x) && any(regexpi(x, '.ns5'))))
    p.addParameter('ns6', '', @(x) (ischar(x) && any(regexpi(x, '.ns6'))))
    p.addParameter('~ch', [], @isnumeric);
    p.addParameter('class', 'double', @ischar)
    p.addParameter('trialdef', {''}, @iscell)

    p.parse(nev, optionalInputs{:});
end