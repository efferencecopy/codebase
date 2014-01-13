function [d,h, wf] = my_abfload(fn)

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
TagInfo=define_TagInfo;



% -------------------------------------------------------------------------
%      read parameters of interest
% -------------------------------------------------------------------------


% read values from header
for g=1:numOfParams
    fseek(fid, headPar(g).offs,'bof');
    sz=length(headPar(g).value);
    h.(headPar(g).name) = fread(fid, sz, headPar(g).numType); % use dynamic field names
end
 
h.fFileSignature = h.fFileSignature'; % transposed
h.fFileVersionNumber = double(h.fFileVersionNumber)' * [0.0001; 0.001; 0.1; 1];
h.lFileStartTime = h.uFileStartTimeMS*.001;  % convert ms to s


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

% this is a hack: determine where either of strings begin
progString='clampex';
goodstart= regexpi(char(BigString)',progString);
if isempty(goodstart); warning('problems in StringsSection'); end
BigString=BigString(goodstart(1):end)';
stringends=find(BigString==0);
stringends=[0 stringends];
for i=1:length(stringends)-1
    Strings{i}=char(BigString(stringends(i)+1:stringends(i+1)-1));
end


h.recChNames=[];
h.recChUnits=[];

% --- read in the ADCSection & copy some values to header h
for i=1:ADCSection.llNumEntries
    ADCsec=ReadSection(fid,ADCSection.uBlockIndex*BLOCKSIZE+ADCSection.uBytes*(i-1),ADCInfo);
    ii=ADCsec.nADCNum+1;
    h.nADCSamplingSeq(i)=ADCsec.nADCNum;
    h.recChNames=strvcat(h.recChNames, Strings{ADCsec.lADCChannelNameIndex});
    unitsIndex=ADCsec.lADCUnitsIndex;
    if unitsIndex>0
        h.recChUnits=strvcat(h.recChUnits, Strings{ADCsec.lADCUnitsIndex});
    else
        h.recChUnits=strvcat(h.recChUnits,'');
    end
    h.nTelegraphEnable(ii)=ADCsec.nTelegraphEnable;
    h.fTelegraphAdditGain(ii)=ADCsec.fTelegraphAdditGain;
    h.fInstrumentScaleFactor(ii)=ADCsec.fInstrumentScaleFactor;
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
h.nOperationMode=ProtocolSec.nOperationMode;
h.fSynchTimeUnit=ProtocolSec.fSynchTimeUnit;
h.nExperimentType = ProtocolSec.nExperimentType;
h.nADCNumChannels=ADCSection.llNumEntries;
h.lActualAcqLength=DataSection.llNumEntries;
h.lDataSectionPtr=DataSection.uBlockIndex;
h.nNumPointsIgnored=0;

% in ABF version < 2.0 h.fADCSampleInterval is the sampling interval
% defined as
%     1/(sampling freq*number_of_channels)
% so divide ProtocolSec.fADCSequenceInterval by the number of channels
h.fADCSampleInterval=ProtocolSec.fADCSequenceInterval/h.nADCNumChannels;
h.fADCRange=ProtocolSec.fADCRange;
h.lADCResolution=ProtocolSec.lADCResolution;

% --- in contrast to procedures with all other sections do not read the
% sync array section but rather copy the values of its fields to the
% corresponding fields of h
lSynchArrayPtr=SynchArraySection.uBlockIndex;
lSynchArraySize=SynchArraySection.llNumEntries;




% -------------------------------------------------------------------------
%     groom parameters & perform some plausibility checks
% -------------------------------------------------------------------------
if h.lActualAcqLength<h.nADCNumChannels,
    fclose(fid);
    error('less data points than sampled channels in file');
end
% the numerical value of all recorded channels (numbers 0..15)
recChIdx=h.nADCSamplingSeq(1:h.nADCNumChannels);
% the corresponding indices into loaded data d
recChInd=1:length(recChIdx);

% convert to cell arrays
h.recChNames=deblank(cellstr(h.recChNames));
h.recChUnits=deblank(cellstr(h.recChUnits));


% gain of telegraphed instruments, if any
addGain=h.nTelegraphEnable.*h.fTelegraphAdditGain;
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
headOffset=h.lDataSectionPtr*BLOCKSIZE+h.nNumPointsIgnored*dataSz;

% h.fADCSampleInterval is the TOTAL samp interval for a single channel (?)
% not the time between samples for a single channel... Define the sample
% interval accordingly
h.si=h.fADCSampleInterval*h.nADCNumChannels;

nSweeps=h.lActualEpisodes;
sweeps=1:h.lActualEpisodes;

% determine time unit in synch array section
switch h.fSynchTimeUnit
    case 0
        % time information in synch array section is in terms of ticks
        h.synchArrTimeBase=1;
    otherwise
        % time information in synch array section is in terms of usec
        h.synchArrTimeBase=h.fSynchTimeUnit;
end

% read in the TagSection, do a few computations & write to h.tags
h.tags=[];
for i=1:TagSection.llNumEntries
    tmp=ReadSection(fid,TagSection.uBlockIndex*BLOCKSIZE+TagSection.uBytes*(i-1),TagInfo);
    % time of tag entry from start of experiment in s (corresponding expisode
    % number, if applicable, will be determined later)
    h.tags(i).timeSinceRecStart=tmp.lTagTime*h.synchArrTimeBase/1e6;
    h.tags(i).comment=char(tmp.sComment)';
end




% -------------------------------------------------------------------------
%     read data (note: from here on code is generic and abf version
%    should not matter)
% -------------------------------------------------------------------------
switch h.nOperationMode
        
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
        h.sweepLengthInPts=synchArr(1,2)/h.nADCNumChannels;
        
        % the starting ticks of episodes in sample points (t0=1=beginning of
        % recording)
        h.sweepStartInPts=synchArr(:,1)*(h.synchArrTimeBase/h.fADCSampleInterval/h.nADCNumChannels);
        
        
        % determine first point and number of points to be read
        startPt=0;
        h.dataPts=h.lActualAcqLength;
        h.dataPtsPerChan=h.dataPts/h.nADCNumChannels;
        if rem(h.dataPts,h.nADCNumChannels)>0 || rem(h.dataPtsPerChan,h.lActualEpisodes)>0
            fclose(fid);
            error('number of data points not OK');
        end
        
        % temporary helper var
        dataPtsPerSweep=h.sweepLengthInPts*h.nADCNumChannels;
        d=zeros(h.sweepLengthInPts,length(recChInd),nSweeps);
        
        % the starting ticks of episodes in sample points WITHIN THE DATA FILE
        selectedSegStartInPts=((sweeps-1)*dataPtsPerSweep)*dataSz+headOffset;
        
        % ** load data
        fseek(fid,startPt*dataSz+headOffset,'bof'); 
        for i=1:nSweeps,
            fseek(fid,selectedSegStartInPts(i),'bof');
            [tmpd,n]=fread(fid,dataPtsPerSweep,precision);
            h.dataPtsPerChan=n/h.nADCNumChannels;
            if rem(n,h.nADCNumChannels)>0
                fclose(fid);
                error('number of data points in episode not OK');
            end
            
            % separate channels..
            tmpd=reshape(tmpd,h.nADCNumChannels,h.dataPtsPerChan);
            tmpd=tmpd';
            
            % if data format is integer, scale appropriately; if it's float, d is fine
            if ~h.nDataFormat
                for j=1:length(recChInd),
                    ch=recChIdx(recChInd(j))+1;
                    tmpd(:,j)=tmpd(:,j)/(h.fInstrumentScaleFactor(ch)*h.fSignalGain(ch)*h.fADCProgrammableGain(ch)*addGain(ch))...
                        *h.fADCRange/h.lADCResolution+h.fInstrumentOffset(ch)-h.fSignalOffset(ch);
                end
            end
            
            % now fill 3d array
            d(:,:,i)=tmpd;
        end
        
        %
        % now load the waveforms (CAH added this section)
        %
        % This part doesn't quite work. For the one example file it looks
        % like the command waveform comes after the analog wave form....
        % It's almost like I'm reading the same data but with a small
        % delay. Maybe the 'headOffset_wf' is wrong...
        %
        wf = nan(size(d));
        headOffset_wf = DACSection.uBlockIndex*BLOCKSIZE+h.nNumPointsIgnored*dataSz;
        selectedSegStartInPts=((sweeps-1)*dataPtsPerSweep)*dataSz+headOffset_wf;
        fseek(fid,startPt*dataSz+headOffset_wf,'bof'); 
        for i=1:nSweeps,
            fseek(fid,selectedSegStartInPts(i),'bof');
            [tmpd,n]=fread(fid,dataPtsPerSweep,precision);
            h.dataPtsPerChan=n/h.nADCNumChannels;
            if rem(n,h.nADCNumChannels)>0
                fclose(fid);
                error('number of data points in episode not OK');
            end
            
            % separate channels..
            tmpd=reshape(tmpd,h.nADCNumChannels,h.dataPtsPerChan);
            tmpd=tmpd';
            
            % if data format is integer, scale appropriately; if it's float, d is fine
            if ~h.nDataFormat
                if i==1
                    warning('DAC section probably not being scaled appropriately')
                end
                for j=1:length(recChInd),
                    ch=recChIdx(recChInd(j))+1;
                    tmpd(:,j)=tmpd(:,j)/(h.fInstrumentScaleFactor(ch)*h.fSignalGain(ch)*h.fADCProgrammableGain(ch)*addGain(ch))...
                        *h.fADCRange/h.lADCResolution+h.fInstrumentOffset(ch)-h.fSignalOffset(ch);
                end
            end
            
            % now fill 3d array
            wf(:,:,i)=tmpd;
        end
        %
        % END CH ADDITIONS
        %
        
        
        
    case 3
        % start at the beginning
        startPt = 0;
        
        % define the stop point (take all the data)
        h.dataPtsPerChan=h.lActualAcqLength/h.nADCNumChannels;
        h.dataPts=h.dataPtsPerChan*h.nADCNumChannels;

        
        % read in some data
        fseek(fid,startPt*dataSz+headOffset,'bof');
        d = fread(fid,h.dataPts,precision);
        
        % separate channels..
        d = reshape(d, h.nADCNumChannels, h.dataPtsPerChan);
        d = d';
        
        
        % if data format is integer, scale appropriately; if it's float, d is fine
        if ~h.nDataFormat
            for j=1:length(recChInd),
                ch=recChIdx(recChInd(j))+1;
                d(:,j)=d(:,j)/(h.fInstrumentScaleFactor(ch)*h.fSignalGain(ch)*h.fADCProgrammableGain(ch)*addGain(ch))...
                    *h.fADCRange/h.lADCResolution+h.fInstrumentOffset(ch)-h.fSignalOffset(ch);
            end
        end
        
    otherwise
        disp('unknown recording mode -- returning empty matrix');
        d=[];
        h.si=[];
end
fclose(fid);

% finally, possibly add information on episode number to tags
if ~isempty(h.tags) && isfield(h,'sweepStartInPts')
    for i=1:numel(h.tags)
        tmp=find(h.tags(i).timeSinceRecStart>=h.sweepStartInPts/1e6*h.si);
        h.tags(i).episodeIndex=tmp(end);
    end
end

end % function


% ########################################################################
%                         LOCAL FUNCTIONS
% ########################################################################

function out=define_header
    headPar={
        'fFileSignature',0,'*char',[-1 -1 -1 -1];
        'fFileVersionNumber',4,'bit8=>int',[-1 -1 -1 -1];
        'uFileInfoSize',8,'uint32',-1;
        'lActualEpisodes',12,'uint32',-1;
        'uFileStartDate',16','uint32',-1;
        'uFileStartTimeMS',20,'uint32',-1;
        'uStopwatchTime',24,'uint32',-1;
        'nFileType',28,'int16',-1;
        'nDataFormat',30,'int16',-1;
        'nSimultaneousScan',32,'int16',-1;
        'nCRCEnable',34,'int16',-1;
        'uFileCRC',36,'uint32',-1;
        'FileGUID',40,'uint32',-1;
        'uCreatorVersion',56,'uint32',-1;
        'uCreatorNameIndex',60,'uint32',-1;
        'uModifierVersion',64,'uint32',-1;
        'uModifierNameIndex',68,'uint32',-1;
        'uProtocolPathIndex',72,'uint32',-1;
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

function TagInfo=define_TagInfo
    TagInfo={
        'lTagTime','int32',1;
        'sComment','char',56;
        'nTagType','int16',1;
        'nVoiceTagNumber_or_AnnotationIndex','int16',1;
        };
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

