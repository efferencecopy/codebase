function paradigmName = paradigmLibrary(paradigmID)
%Each rex program should have a unique paradigmID. This code is droped into
%the datafile by rex and is essential for nex2stro to function correctly.
%Nex2stro must have access to the proper codes (specific to each paradigm)
%and the paradigm id maps the rex data file to the appropriate set of
%codes.

switch paradigmID
    case 100
        paradigmName = 'WhiteNoiseCodes';
    case 101
        paradigmName = 'GaborEdgeCodes';
    case 102
        paradigmName = 'SMurrayCodes';
    case 103
        paradigmName = 'NeuroThreshCodes';
    case 104
        paradigmName = 'SMurray1Codes';
    case 105
        paradigmName = 'FixStimCodes';
    case 106
        paradigmName = 'GridLMPlaneCodes';
    case 107
        paradigmName = 'IsoSampCodes';
    case 150
        paradigmName = 'GratingCodes';
    case 200
        paradigmName = 'DTCodesTraining';
    case 210
        paradigmName = 'DTCodes';
    case 211
        paradigmName = 'HABITCodes';
    case 212
        paradigmName = 'DTStimCodes';
    case 213
        paradigmName = 'DTScotCodes';
    case 250
        paradigmName = 'NM2SCodes';
    otherwise
        error('Unable to map paradigmID <%d> to a paradigm header file', paradigmID);
end