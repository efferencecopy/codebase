function [d,h] = my_abfload(fn)

% **  [data, header] = abfload(filename)
%
% >>> INPUT VARIABLES >>>
% NAME        TYPE, DEFAULT      DESCRIPTION
% fn          char array         abf data file name
%
% << OUTPUT VARIABLES <<<
% NAME  TYPE            DESCRIPTION
% d                     the data read, the format depending on the record-
%                        ing mode
%   1. GAP-FREE:
%   2d array        2d array of size <data pts> by <number of chans>
%                    Examples of access:
%                    d(:,2)       data from channel 2 at full length
%                    d(1:100,:)   first 100 data points from all channels
%
%   2. EPISODIC FIXED-LENGTH/WAVEFORM FIXED-LENGTH/HIGH-SPEED OSCILLOSCOPE:
%   3d array        3d array of size <data pts per sweep> by <number of chans> by <number of sweeps>.
%                    Examples of access:
%                    d(:,2,:)            a matrix containing all episodes
%                                        (at full length) of the second
%                                        channel in its columns
%                    d(1:200,:,[1 11])   contains first 200 data points of
%                                        episodes 1 and 11 of all channels
% 
%
% CONTRIBUTORS
%   Original version by Harald Hentschke (harald.hentschke@uni-tuebingen.de)
%   Extended to abf version 2.0 by Forrest Collman (fcollman@Princeton.edu)
%



% define constants
BLOCKSIZE=512;


% -------------------------------------------------------------------------
%                       determine abf version
% -------------------------------------------------------------------------
[fid, messg] = fopen(fn,'r','ieee-le');
if fid == -1, error(messg); end


% determine absolute file size
fseek(fid,0,'eof');
fileSz=ftell(fid);



% *** read value of parameter 'fFileSignature' (i.e. abf version) from header ***
fseek(fid,0,'bof');
fFileSignature = fread(fid, 4,'uint8=>char')';
fseek(fid,0,'bof'); % rewind


% one of the first checks must be whether file signature is valid
if isempty(regexp('ABF2', fFileSignature, 'once'))
    fclose(fid);
    error('unknown or incompatible file signature. ABF file must be 2.0 or newer');
end

% -------------------------------------------------------------------------
%    define file information ('header' parameters) of interest
% -------------------------------------------------------------------------
% The list of header parameters created below (variable 'headPar') is
% derived from the abf version 1.8 header section. It is by no means
% exhaustive (i.e. there are many more parameters in abf files) but
% sufficient for proper upload, scaling and arrangement of data acquired
% under many conditons. Further below, these parameters will be made fields
% of struct h. h, which is also an output variable, is then used in PART 3,
% which does the actual job of uploading, scaling and rearranging the data.
% That part of the code relies on h having a certain set of fields
% irrespective of ABF version.
% Unfortunately, in the transition to ABF version 2.0 many of the header
% parameters were moved to different places within the abf file and/or
% given other names or completely restructured. In order for the code to
% work with pre- and post-2.0 data files, all parameters missing in the
% header must be gotten into h. This is accomplished in lines ~288 and
% following:
%     if h.fFileVersionNumber>=2
%       ...
% Furthermore,
% - h as an output from an ABF version < 2.0 file will not contain new
%   parameters introduced into the header like 'nCRCEnable'
% - h will in any case contain a few 'home-made' fields that have
%   proven to be useful. Some of them depend on the recording mode. Among
%   the more or less self-explanatory ones are
% -- si                   sampling interval
% -- recChNames           the names of all channels, e.g. 'IN 8',...
% -- dataPtsPerChan       sample points per channel
% -- dataPts              sample points in file
% -- recTime              recording start and stop time in seconds from
%                         midnight (millisecond resolution)
% -- sweepLengthInPts     sample points per sweep (one channel)
% -- sweepStartInPts      the start times of sweeps in sample points
%                         (from beginning of recording)


% define header 
headPar=define_header;
numOfParams=size(headPar,1);

% define some other sections sections
Sections=define_Sections;
ProtocolInfo=define_ProtocolInfo;
ADCInfo=define_ADCInfo;


% read values from header
for g=1:numOfParams
    fseek(fid, headPar(g).offs,'bof');
    sz=length(headPar(g).value);
    h.(headPar(g).name) = fread(fid, sz, headPar(g).numType); % use dynamic field names
end
 

% -----------------------------------------------------------------------
% *** read file information that has moved from the header section to
% other sections in ABF version >= 2.0 and assign selected values to
% fields of 'generic' header variable h ***
% -----------------------------------------------------------------------
numOfSections=length(Sections);
offset=76;

% this creates all sections (ADCSection, ProtocolSection, etc.)
for i=1:numOfSections
    eval([Sections(i).name '=ReadSectionInfo(fid,offset);']);
    offset=offset+4+4+8;
end

% --- read in the StringsSection and use some fields (to retrieve
% information on the names of recorded channels and the units). The path to
% the clampex protocol is in Strings{2} later on...
fseek(fid,StringsSection.uBlockIndex*BLOCKSIZE,'bof');
BigString=fread(fid,StringsSection.uBytes,'char');

BigString = char(BigString)';
goodstart= regexpi(BigString, 'clampex', 'once');
assert(~isempty(goodstart), 'LOADABF ERROR: problem in strings section');

BigString=BigString(goodstart:end);
stringends=regexpi(BigString, char(0));
stringends=[0 stringends];
for i=1:length(stringends)-1
    Strings{i}=BigString(stringends(i)+1:stringends(i+1)-1);
end

% retrieve the protocol name
fseps = regexpi(Strings{2}, filesep);
h.protocolName = Strings{2}(fseps(end)+1:end);



% --- read in the ADCSection & copy some values to header h
h.recChNames=[];
h.recChUnits=[];
for i=1:ADCSection.llNumEntries
    ADCsec=ReadSection(fid,ADCSection.uBlockIndex*BLOCKSIZE+ADCSection.uBytes*(i-1),ADCInfo);
    ii=ADCsec.nADCNum+1;
    nADCSamplingSeq(i)=ADCsec.nADCNum;
    h.recChNames=strvcat(h.recChNames, Strings{ADCsec.lADCChannelNameIndex});
    unitsIndex=ADCsec.lADCUnitsIndex;
    if unitsIndex>0
        h.recChUnits=strvcat(h.recChUnits, Strings{ADCsec.lADCUnitsIndex});
    else
        h.recChUnits=strvcat(h.recChUnits,'');
    end
    nTelegraphEnable(ii)=ADCsec.nTelegraphEnable;
    fTelegraphAdditGain(ii)=ADCsec.fTelegraphAdditGain;
    fInstrumentScaleFactor(ii)=ADCsec.fInstrumentScaleFactor;
    h.fSignalGain(ii)=ADCsec.fSignalGain;
    h.fADCProgrammableGain(ii)=ADCsec.fADCProgrammableGain;
    h.fInstrumentOffset(ii)=ADCsec.fInstrumentOffset;
    h.fSignalOffset(ii)=ADCsec.fSignalOffset;
    h.fTelegraphFilter(ii) = ADCsec.fTelegraphFilter;
    h.fSignalLowpassFilter(ii) = ADCsec.fSignalLowpassFilter;
    h.fSignalHighpassFilter(ii) = ADCsec.fSignalHighpassFilter;
end


% --- read in the protocol section & copy some values to header h
ProtocolSec=ReadSection(fid,ProtocolSection.uBlockIndex*BLOCKSIZE,ProtocolInfo);
fSynchTimeUnit=ProtocolSec.fSynchTimeUnit;
h.nExperimentType = ProtocolSec.nExperimentType;
nADCNumChannels=ADCSection.llNumEntries;
lActualAcqLength=DataSection.llNumEntries;
nNumPointsIgnored=0;

% in ABF version < 2.0 fADCSampleInterval is the sampling interval
% defined as
%     1/(sampling freq*number_of_channels)
% so divide ProtocolSec.fADCSequenceInterval by the number of channels
fADCSampleInterval=ProtocolSec.fADCSequenceInterval/nADCNumChannels;

% --- in contrast to procedures with all other sections do not read the
% sync array section but rather copy the values of its fields to the
% corresponding fields of h
lSynchArrayPtr=SynchArraySection.uBlockIndex;
lSynchArraySize=SynchArraySection.llNumEntries;




% -------------------------------------------------------------------------
%     groom parameters & perform some plausibility checks
% -------------------------------------------------------------------------
if lActualAcqLength<nADCNumChannels,
    fclose(fid);
    error('less data points than sampled channels in file');
end
% the numerical value of all recorded channels (numbers 0..15)
recChIdx=nADCSamplingSeq(1:nADCNumChannels);
% the corresponding indices into loaded data d
recChInd=1:length(recChIdx);

% convert to cell arrays
h.recChNames=deblank(cellstr(h.recChNames));
h.recChUnits=deblank(cellstr(h.recChUnits));


% gain of telegraphed instruments, if any
addGain=nTelegraphEnable.*fTelegraphAdditGain;
addGain(addGain==0)=1;


% determine offset at which data start
switch h.nDataFormat
    case 0
        dataSz=2;  % bytes/point
        precision='int16';
    case 1
        dataSz=4;  % bytes/point
        precision='float32';
    otherwise
        fclose(fid);
        error('invalid number format');
end
headOffset=DataSection.uBlockIndex*BLOCKSIZE+nNumPointsIgnored*dataSz;

% fADCSampleInterval is the TOTAL samp interval for a single channel (?)
% not the time between samples for a single channel... Define the sample
% interval accordingly
sampInt=fADCSampleInterval*nADCNumChannels;

nSweeps=h.lActualEpisodes;
sweeps=1:nSweeps;

% determine time unit in synch array section
switch fSynchTimeUnit
    case 0
        % time information in synch array section is in terms of ticks
        synchArrTimeBase=1;
    otherwise
        % time information in synch array section is in terms of usec
        synchArrTimeBase=fSynchTimeUnit;
end



switch ProtocolSec.nOperationMode
        
    case {2,4,5}  % 2=> event-driven fixed-length; 4=> high speed osciliscope 5=> waveform fixed-length

        % extract timing information on sweeps
        if (lSynchArrayPtr<=0 || lSynchArraySize<=0),
            fclose(fid);
            error('internal variables ''lSynchArraynnn'' are zero or negative');
        end
        % the byte offset at which the SynchArraySection starts
        lSynchArrayPtrByte=BLOCKSIZE*lSynchArrayPtr;
        
        % before reading Synch Arr parameters check if file is big enough to hold them
        % 4 bytes/long, 2 values per episode (start and length)
        if lSynchArrayPtrByte+2*4*lSynchArraySize>fileSz,
            fclose(fid);
            error('file seems not to contain complete Synch Array Section');
        end
        
        
        fseek(fid,lSynchArrayPtrByte,'bof');
        synchArr=fread(fid,lSynchArraySize*2,'int32');

        % make synchArr a lSynchArraySize x 2 matrix
        synchArr=permute(reshape(synchArr',2,lSynchArraySize),[2 1]);
        if numel(unique(synchArr(:,2)))>1
            fclose(fid);
            error('sweeps of unequal length in file recorded in fixed-length mode');
        end
        
        % the length of sweeps in sample points (**note: parameter lLength of
        % the ABF synch section is expressed in samples (ticks) whereas
        % parameter lStart is given in synchArrTimeBase units)
        sweepLengthInPts=synchArr(1,2)/nADCNumChannels;
        
        % the starting ticks of episodes in sample points (t0=1=beginning of
        % recording)
        sweepStartInPts=synchArr(:,1)*(synchArrTimeBase/fADCSampleInterval/nADCNumChannels);
        
        
        % determine first point and number of points to be read
        startPt=0;
        nDataPts=lActualAcqLength;
        nDataPtsPerChan=nDataPts/nADCNumChannels;
        if rem(nDataPts,nADCNumChannels)>0 || rem(nDataPtsPerChan,h.lActualEpisodes)>0
            fclose(fid);
            error('number of data points not OK');
        end
        
        % temporary helper var
        dataPtsPerSweep=sweepLengthInPts*nADCNumChannels;
        d=zeros(sweepLengthInPts,length(recChInd),nSweeps);
        
        % the starting ticks of episodes in sample points WITHIN THE DATA FILE
        selectedSegStartInPts=((sweeps-1)*dataPtsPerSweep)*dataSz+headOffset;
        
        % ** load data
        fseek(fid,startPt*dataSz+headOffset,'bof'); 
        for i=1:nSweeps,
            fseek(fid,selectedSegStartInPts(i),'bof');
            [tmpd,n]=fread(fid,dataPtsPerSweep,precision);
            nDataPtsPerChan=n/nADCNumChannels;
            if rem(n,nADCNumChannels)>0
                fclose(fid);
                error('number of data points in episode not OK');
            end
            
            % separate channels..
            tmpd=reshape(tmpd,nADCNumChannels,nDataPtsPerChan);
            tmpd=tmpd';
            
            % if data format is integer, scale appropriately; if it's float, d is fine
            if ~h.nDataFormat
                for j=1:length(recChInd),
                    ch=recChIdx(recChInd(j))+1;
                    tmpd(:,j)=tmpd(:,j)/(fInstrumentScaleFactor(ch)*h.fSignalGain(ch)*h.fADCProgrammableGain(ch)*addGain(ch))...
                        *ProtocolSec.fADCRange/ProtocolSec.lADCResolution+h.fInstrumentOffset(ch)-h.fSignalOffset(ch);
                end
            end
            
            % now fill 3d array
            d(:,:,i)=tmpd;
        end
        
               
    case 3
        % start at the beginning
        startPt = 0;
        
        % define the stop point (take all the data)
        nDataPtsPerChan=lActualAcqLength/nADCNumChannels;
        nDataPts=nDataPtsPerChan*nADCNumChannels;

        
        % read in some data
        fseek(fid,startPt*dataSz+headOffset,'bof');
        d = fread(fid,nDataPts,precision);
        
        % separate channels..
        d = reshape(d, nADCNumChannels, nDataPtsPerChan);
        d = d';
        
        
        % if data format is integer, scale appropriately; if it's float, d is fine
        if ~h.nDataFormat
            for j=1:length(recChInd),
                ch=recChIdx(recChInd(j))+1;
                d(:,j)=d(:,j)/(fInstrumentScaleFactor(ch)*h.fSignalGain(ch)*h.fADCProgrammableGain(ch)*addGain(ch))...
                    *ProtocolSec.fADCRange/ProtocolSec.lADCResolution+h.fInstrumentOffset(ch)-h.fSignalOffset(ch);
            end
        end
        
    otherwise
        disp('unknown recording mode -- returning empty matrix');
        d=[];
        h=[];
end

% technically I could just close fid, but closing 'all' avoids there being
% open files unbeknownst to the user, which could interfere with the
% ability to open them simultaneously in matlab and clampfit.
fclose('all');

end % function


% ########################################################################
%                         LOCAL FUNCTIONS
% ########################################################################

function out=define_header
    
    % NOTE: C.Hass commented out most of the options below because I don't
    % know what they do, they don't seem to be productive, and I find them
    % annoying (2/20/2014)
    
    headPar={
        %'fFileSignature',0,'*char',[-1 -1 -1 -1];
        %'fFileVersionNumber',4,'bit8=>int',[-1 -1 -1 -1];
        %'uFileInfoSize',8,'uint32',-1;
        'lActualEpisodes',12,'uint32',-1;
        %'uFileStartDate',16','uint32',-1;
        %'uFileStartTimeMS',20,'uint32',-1;
        %'uStopwatchTime',24,'uint32',-1;
        %'nFileType',28,'int16',-1;
        'nDataFormat',30,'int16',-1;
        %'nSimultaneousScan',32,'int16',-1;
        %'nCRCEnable',34,'int16',-1;
        %'uFileCRC',36,'uint32',-1;
        %'FileGUID',40,'uint32',-1;
        %'uCreatorVersion',56,'uint32',-1;
        %'uCreatorNameIndex',60,'uint32',-1;
        %'uModifierVersion',64,'uint32',-1;
        %'uModifierNameIndex',68,'uint32',-1;
        %'uProtocolPathIndex',72,'uint32',-1;
        };
    out = cell2struct(headPar,{'name','offs','numType','value'},2);
end

function out=define_Sections
    Sections={'ProtocolSection';
        'ADCSection';
        'DACSection';
        'EpochSection';
        'ADCPerDACSection';
        'EpochPerDACSection';
        'UserListSection';
        'StatsRegionSection';
        'MathSection';
        'StringsSection';
        'DataSection';
        'TagSection';
        'ScopeSection';
        'DeltaSection';
        'VoiceTagSection';
        'SynchArraySection';
        'AnnotationSection';
        'StatsSection';
        };
    out = cell2struct(Sections,{'name'},2);
end

function ProtocolInfo=define_ProtocolInfo
    ProtocolInfo={
        'nOperationMode','int16',1;
        'fADCSequenceInterval','float',1;
        'bEnableFileCompression','bit1',1;
        'sUnused1','char',3;
        'uFileCompressionRatio','uint32',1;
        'fSynchTimeUnit','float',1;
        'fSecondsPerRun','float',1;
        'lNumSamplesPerEpisode','int32',1;
        'lPreTriggerSamples','int32',1;
        'lEpisodesPerRun','int32',1;
        'lRunsPerTrial','int32',1;
        'lNumberOfTrials','int32',1;
        'nAveragingMode','int16',1;
        'nUndoRunCount','int16',1;
        'nFirstEpisodeInRun','int16',1;
        'fTriggerThreshold','float',1;
        'nTriggerSource','int16',1;
        'nTriggerAction','int16',1;
        'nTriggerPolarity','int16',1;
        'fScopeOutputInterval','float',1;
        'fEpisodeStartToStart','float',1;
        'fRunStartToStart','float',1;
        'lAverageCount','int32',1;
        'fTrialStartToStart','float',1;
        'nAutoTriggerStrategy','int16',1;
        'fFirstRunDelayS','float',1;
        'nChannelStatsStrategy','int16',1;
        'lSamplesPerTrace','int32',1;
        'lStartDisplayNum','int32',1;
        'lFinishDisplayNum','int32',1;
        'nShowPNRawData','int16',1;
        'fStatisticsPeriod','float',1;
        'lStatisticsMeasurements','int32',1;
        'nStatisticsSaveStrategy','int16',1;
        'fADCRange','float',1;
        'fDACRange','float',1;
        'lADCResolution','int32',1;
        'lDACResolution','int32',1;
        'nExperimentType','int16',1;
        'nManualInfoStrategy','int16',1;
        'nCommentsEnable','int16',1;
        'lFileCommentIndex','int32',1;
        'nAutoAnalyseEnable','int16',1;
        'nSignalType','int16',1;
        'nDigitalEnable','int16',1;
        'nActiveDACChannel','int16',1;
        'nDigitalHolding','int16',1;
        'nDigitalInterEpisode','int16',1;
        'nDigitalDACChannel','int16',1;
        'nDigitalTrainActiveLogic','int16',1;
        'nStatsEnable','int16',1;
        'nStatisticsClearStrategy','int16',1;
        'nLevelHysteresis','int16',1;
        'lTimeHysteresis','int32',1;
        'nAllowExternalTags','int16',1;
        'nAverageAlgorithm','int16',1;
        'fAverageWeighting','float',1;
        'nUndoPromptStrategy','int16',1;
        'nTrialTriggerSource','int16',1;
        'nStatisticsDisplayStrategy','int16',1;
        'nExternalTagType','int16',1;
        'nScopeTriggerOut','int16',1;
        'nLTPType','int16',1;
        'nAlternateDACOutputState','int16',1;
        'nAlternateDigitalOutputState','int16',1;
        'fCellID','float',3;
        'nDigitizerADCs','int16',1;
        'nDigitizerDACs','int16',1;
        'nDigitizerTotalDigitalOuts','int16',1;
        'nDigitizerSynchDigitalOuts','int16',1;
        'nDigitizerType','int16',1;
        };
end

function ADCInfo=define_ADCInfo
    ADCInfo={
        'nADCNum','int16',1;
        'nTelegraphEnable','int16',1;
        'nTelegraphInstrument','int16',1;
        'fTelegraphAdditGain','float',1;
        'fTelegraphFilter','float',1;
        'fTelegraphMembraneCap','float',1;
        'nTelegraphMode','int16',1;
        'fTelegraphAccessResistance','float',1;
        'nADCPtoLChannelMap','int16',1;
        'nADCSamplingSeq','int16',1;
        'fADCProgrammableGain','float',1;
        'fADCDisplayAmplification','float',1;
        'fADCDisplayOffset','float',1;
        'fInstrumentScaleFactor','float',1;
        'fInstrumentOffset','float',1;
        'fSignalGain','float',1;
        'fSignalOffset','float',1;
        'fSignalLowpassFilter','float',1;
        'fSignalHighpassFilter','float',1;
        'nLowpassFilterType','char',1;
        'nHighpassFilterType','char',1;
        'fPostProcessLowpassFilter','float',1;
        'nPostProcessLowpassFilterType','char',1;
        'bEnabledDuringPN','bit1',1;
        'nStatsChannelPolarity','int16',1;
        'lADCChannelNameIndex','int32',1;
        'lADCUnitsIndex','int32',1;
        };
end

function EpochInfo = define_EpochInfo

EpochInfoPerDACDescription = {
       'nEpochNum', 'int16', 1;
       'nDACNum', 'int16', 1;
       'nEpochType', 'int16', 1;
       'fEpochInitLevel' ,'float', 1;
       'fEpochLevelInc', 'float', 1;
       'lEpochInitDuration', 'uint32', 1;
       'lEpochDurationInc', 'uint32', 1;
       'lEpochPulsePeriod', 'uint32', 1;
       'lEpochPulseWidth', 'uint32', 1;
       'sUnused', 'char', 18};

EpochInfoDescription = {
       'nEpochNum','int16', 1;
       'nDigitalValue','int16', 1;
       'nDigitalTrainValue','int16', 1;
       'nAlternateDigitalValue','int16', 1;
       'nAlternateDigitalTrainValue','int16', 1;
       'bEpochCompression','bit1', 1;
       'sUnused','char', 21};
end

function Section=ReadSection(fid,offset,Format)
    s=cell2struct(Format,{'name','numType','number'},2);
    fseek(fid,offset,'bof');
    for i=1:length(s)
        eval(['[Section.' s(i).name ',n]=fread(fid,' num2str(s(i).number) ',''' s(i).numType ''');']);
    end
end

function SectionInfo=ReadSectionInfo(fid,offset)
    fseek(fid,offset,'bof');
    SectionInfo.uBlockIndex=fread(fid,1,'uint32');
    fseek(fid,offset+4,'bof');
    SectionInfo.uBytes=fread(fid,1,'uint32');
    fseek(fid,offset+8,'bof');
    SectionInfo.llNumEntries=fread(fid,1,'int64');
end


%#### SOME PYTHON CODE I GRABBED ONLINE. EPOCH RELATED THINGS MIGHT CONTAIN
%#### WAVEFORMS, BUT I'LL NEED TO USE OTHER SECTIONS AS A ROSETTA STONE TO
%#### FIURE OUT DATA TYPES...
%#### More details here: https://github.com/NeuralEnsemble/python-neo/blob/master/neo/io/axonio.py

% % headerDescriptionV1= [
% %          ('fFileSignature',0,'4s'),
% %          ('fFileVersionNumber',4,'f' ),
% %          ('nOperationMode',8,'h' ),
% %          ('lActualAcqLength',10,'i' ),
% %          ('nNumPointsIgnored',14,'h' ),
% %          ('lActualEpisodes',16,'i' ),
% %          ('lFileStartTime',24,'i' ),
% %          ('lDataSectionPtr',40,'i' ),
% %          ('lTagSectionPtr',44,'i' ),
% %          ('lNumTagEntries',48,'i' ),
% %          ('lSynchArrayPtr',92,'i' ),
% %          ('lSynchArraySize',96,'i' ),
% %          ('nDataFormat',100,'h' ),
% %          ('nADCNumChannels', 120, 'h'),
% %          ('fADCSampleInterval',122,'f'),
% %          ('fSynchTimeUnit',130,'f' ),
% %          ('lNumSamplesPerEpisode',138,'i' ),
% %          ('lPreTriggerSamples',142,'i' ),
% %          ('lEpisodesPerRun',146,'i' ),
% %          ('fADCRange', 244, 'f' ),
% %          ('lADCResolution', 252, 'i'),
% %          ('nFileStartMillisecs', 366, 'h'),
% %          ('nADCPtoLChannelMap', 378, '16h'),
% %          ('nADCSamplingSeq', 410, '16h'),
% %          ('sADCChannelName',442, '10s'*16),
% %          ('sADCUnits',602, '8s'*16) ,
% %          ('fADCProgrammableGain', 730, '16f'),
% %          ('fInstrumentScaleFactor', 922, '16f'),
% %          ('fInstrumentOffset', 986, '16f'),
% %          ('fSignalGain', 1050, '16f'),
% %          ('fSignalOffset', 1114, '16f'),
% %          ('nTelegraphEnable',4512, '16h'),
% %          ('fTelegraphAdditGain',4576,'16f'),
% %          ]
% % 
% % 
% % headerDescriptionV2 =[
% %          ('fFileSignature',0,'4s' ),
% %          ('fFileVersionNumber',4,'4b') ,
% %          ('uFileInfoSize',8,'I' ) ,
% %          ('lActualEpisodes',12,'I' ) ,
% %          ('uFileStartDate',16,'I' ) ,
% %          ('uFileStartTimeMS',20,'I' ) ,
% %          ('uStopwatchTime',24,'I' ) ,
% %          ('nFileType',28,'H' ) ,
% %          ('nDataFormat',30,'H' ) ,
% %          ('nSimultaneousScan',32,'H' ) ,
% %          ('nCRCEnable',34,'H' ) ,
% %          ('uFileCRC',36,'I' ) ,
% %          ('FileGUID',40,'I' ) ,
% %          ('uCreatorVersion',56,'I' ) ,
% %          ('uCreatorNameIndex',60,'I' ) ,
% %          ('uModifierVersion',64,'I' ) ,
% %          ('uModifierNameIndex',68,'I' ) ,
% %          ('uProtocolPathIndex',72,'I' ) ,
% %          ]
% % 
% % 
% % sectionNames= ['ProtocolSection',
% %              'ADCSection',
% %              'DACSection',
% %              'EpochSection',
% %              'ADCPerDACSection',
% %              'EpochPerDACSection',
% %              'UserListSection',
% %              'StatsRegionSection',
% %              'MathSection',
% %              'StringsSection',
% %              'DataSection',
% %              'TagSection',
% %              'ScopeSection',
% %              'DeltaSection',
% %              'VoiceTagSection',
% %              'SynchArraySection',
% %              'AnnotationSection',
% %              'StatsSection',
% %              ]
% % 
% % 
% % protocolInfoDescription = [
% %          ('nOperationMode','h'),
% %          ('fADCSequenceInterval','f'),
% %          ('bEnableFileCompression','b'),
% %          ('sUnused1','3s'),
% %          ('uFileCompressionRatio','I'),
% %          ('fSynchTimeUnit','f'),
% %          ('fSecondsPerRun','f'),
% %          ('lNumSamplesPerEpisode','i'),
% %          ('lPreTriggerSamples','i'),
% %          ('lEpisodesPerRun','i'),
% %          ('lRunsPerTrial','i'),
% %          ('lNumberOfTrials','i'),
% %          ('nAveragingMode','h'),
% %          ('nUndoRunCount','h'),
% %          ('nFirstEpisodeInRun','h'),
% %          ('fTriggerThreshold','f'),
% %          ('nTriggerSource','h'),
% %          ('nTriggerAction','h'),
% %          ('nTriggerPolarity','h'),
% %          ('fScopeOutputInterval','f'),
% %          ('fEpisodeStartToStart','f'),
% %          ('fRunStartToStart','f'),
% %          ('lAverageCount','i'),
% %          ('fTrialStartToStart','f'),
% %          ('nAutoTriggerStrategy','h'),
% %          ('fFirstRunDelayS','f'),
% %          ('nChannelStatsStrategy','h'),
% %          ('lSamplesPerTrace','i'),
% %          ('lStartDisplayNum','i'),
% %          ('lFinishDisplayNum','i'),
% %          ('nShowPNRawData','h'),
% %          ('fStatisticsPeriod','f'),
% %          ('lStatisticsMeasurements','i'),
% %          ('nStatisticsSaveStrategy','h'),
% %          ('fADCRange','f'),
% %          ('fDACRange','f'),
% %          ('lADCResolution','i'),
% %          ('lDACResolution','i'),
% %          ('nExperimentType','h'),
% %          ('nManualInfoStrategy','h'),
% %          ('nCommentsEnable','h'),
% %          ('lFileCommentIndex','i'),
% %          ('nAutoAnalyseEnable','h'),
% %          ('nSignalType','h'),
% %          ('nDigitalEnable','h'),
% %          ('nActiveDACChannel','h'),
% %          ('nDigitalHolding','h'),
% %          ('nDigitalInterEpisode','h'),
% %          ('nDigitalDACChannel','h'),
% %          ('nDigitalTrainActiveLogic','h'),
% %          ('nStatsEnable','h'),
% %          ('nStatisticsClearStrategy','h'),
% %          ('nLevelHysteresis','h'),
% %          ('lTimeHysteresis','i'),
% %          ('nAllowExternalTags','h'),
% %          ('nAverageAlgorithm','h'),
% %          ('fAverageWeighting','f'),
% %          ('nUndoPromptStrategy','h'),
% %          ('nTrialTriggerSource','h'),
% %          ('nStatisticsDisplayStrategy','h'),
% %          ('nExternalTagType','h'),
% %          ('nScopeTriggerOut','h'),
% %          ('nLTPType','h'),
% %          ('nAlternateDACOutputState','h'),
% %          ('nAlternateDigitalOutputState','h'),
% %          ('fCellID','3f'),
% %          ('nDigitizerADCs','h'),
% %          ('nDigitizerDACs','h'),
% %          ('nDigitizerTotalDigitalOuts','h'),
% %          ('nDigitizerSynchDigitalOuts','h'),
% %          ('nDigitizerType','h'),
% %          ]
% % 
% % 
% % ADCInfoDescription = [
% %          ('nADCNum','h'),
% %          ('nTelegraphEnable','h'),
% %          ('nTelegraphInstrument','h'),
% %          ('fTelegraphAdditGain','f'),
% %          ('fTelegraphFilter','f'),
% %          ('fTelegraphMembraneCap','f'),
% %          ('nTelegraphMode','h'),
% %          ('fTelegraphAccessResistance','f'),
% %          ('nADCPtoLChannelMap','h'),
% %          ('nADCSamplingSeq','h'),
% %          ('fADCProgrammableGain','f'),
% %          ('fADCDisplayAmplification','f'),
% %          ('fADCDisplayOffset','f'),
% %          ('fInstrumentScaleFactor','f'),
% %          ('fInstrumentOffset','f'),
% %          ('fSignalGain','f'),
% %          ('fSignalOffset','f'),
% %          ('fSignalLowpassFilter','f'),
% %          ('fSignalHighpassFilter','f'),
% %          ('nLowpassFilterType','b'),
% %          ('nHighpassFilterType','b'),
% %          ('fPostProcessLowpassFilter','f'),
% %          ('nPostProcessLowpassFilterType','c'),
% %          ('bEnabledDuringPN','b'),
% %          ('nStatsChannelPolarity','h'),
% %          ('lADCChannelNameIndex','i'),
% %          ('lADCUnitsIndex','i'),
% %          ]
% % 
% % TagInfoDescription = [
% %        ('lTagTime','i'),
% %        ('sComment','56s'),
% %        ('nTagType','h'),
% %        ('nVoiceTagNumber_or_AnnotationIndex','h'),
% %        ]
% % 
% % DACInfoDescription = [
% %        ('nDACNum','h'),
% %        ('nTelegraphDACScaleFactorEnable','h'),
% %        ('fInstrumentHoldingLevel', 'f'),
% %        ('fDACScaleFactor','f'),
% %        ('fDACHoldingLevel','f'),
% %        ('fDACCalibrationFactor','f'),
% %        ('fDACCalibrationOffset','f'),
% %        ('lDACChannelNameIndex','i'),
% %        ('lDACChannelUnitsIndex','i'),
% %        ('lDACFilePtr','i'),
% %        ('lDACFileNumEpisodes','i'),
% %        ('nWaveformEnable','h'),
% %        ('nWaveformSource','h'),
% %        ('nInterEpisodeLevel','h'),
% %        ('fDACFileScale','f'),
% %        ('fDACFileOffset','f'),
% %        ('lDACFileEpisodeNum','i'),
% %        ('nDACFileADCNum','h'),
% %        ('nConditEnable','h'),
% %        ('lConditNumPulses','i'),
% %        ('fBaselineDuration','f'),
% %        ('fBaselineLevel','f'),
% %        ('fStepDuration','f'),
% %        ('fStepLevel','f'),
% %        ('fPostTrainPeriod','f'),
% %        ('fPostTrainLevel','f'),
% %        ('nMembTestEnable','h'),
% %        ('nLeakSubtractType','h'),
% %        ('nPNPolarity','h'),
% %        ('fPNHoldingLevel','f'),
% %        ('nPNNumADCChannels','h'),
% %        ('nPNPosition','h'),
% %        ('nPNNumPulses','h'),
% %        ('fPNSettlingTime','f'),
% %        ('fPNInterpulse','f'),
% %        ('nLTPUsageOfDAC','h'),
% %        ('nLTPPresynapticPulses','h'),
% %        ('lDACFilePathIndex','i'),
% %        ('fMembTestPreSettlingTimeMS','f'),
% %        ('fMembTestPostSettlingTimeMS','f'),
% %        ('nLeakSubtractADCIndex','h'),
% %        ('sUnused','124s'),
% %        ]
% % 
% % EpochInfoPerDACDescription = [
% %        ('nEpochNum','h'),
% %        ('nDACNum','h'),
% %        ('nEpochType','h'),
% %        ('fEpochInitLevel','f'),
% %        ('fEpochLevelInc','f'),
% %        ('lEpochInitDuration','i'),
% %        ('lEpochDurationInc','i'),
% %        ('lEpochPulsePeriod','i'),
% %        ('lEpochPulseWidth','i'),
% %        ('sUnused','18s'),
% %        ]
% % 
% % EpochInfoDescription = [
% %        ('nEpochNum','h'),
% %        ('nDigitalValue','h'),
% %        ('nDigitalTrainValue','h'),
% %        ('nAlternateDigitalValue','h'),
% %        ('nAlternateDigitalTrainValue','h'),
% %        ('bEpochCompression','b'),
% %        ('sUnused','21s'),
% %        ]

% % % #### from abf header:
% % 
% % //
% % // Constants for nEpochType
% % //
% % #define ABF_EPOCHDISABLED           0     // disabled epoch
% % #define ABF_EPOCHSTEPPED            1     // stepped waveform
% % #define ABF_EPOCHRAMPED             2     // ramp waveform
% % #define ABF_EPOCH_TYPE_RECTANGLE    3     // rectangular pulse train
% % #define ABF_EPOCH_TYPE_TRIANGLE     4     // triangular waveform
% % #define ABF_EPOCH_TYPE_COSINE       5     // cosinusoidal waveform
% % #define ABF_EPOCH_TYPE_UNUSED       6     // was ABF_EPOCH_TYPE_RESISTANCE
% % #define ABF_EPOCH_TYPE_BIPHASIC     7     // biphasic pulse train
% % #define ABF_EPOCHSLOPE              8     // IonWorks style ramp waveform
% % 
% % % code to extract waveform information:
% % for i=1:EpochPerDACSection.llNumEntries
% %     epoch=ReadSection(fid,EpochPerDACSection.uBlockIndex*BLOCKSIZE+EpochPerDACSection.uBytes*(i-1),EpochInfoPerDACDescription);
% %     epoch
% % end

