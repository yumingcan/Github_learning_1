function [data,si,fhead]=abfload(fn,varargin)
% ** function [d,si,h]=abfload(fn,varargin)
% loads and returns data in ABF (Axon Binary File) format.
% Data may have been acquired in the following modes:
% (1) event-driven variable-length (currently only abf versions < 2.0)
% (2) event-driven fixed-length or waveform-fixed length
% (3) gap-free
% Information about scaling, the time base and the number of channels and 
% episodes is extracted from the header of the abf file.
%
% OPERATION
% If the second input variable is the char array 'info' as in 
%         [d,si,h]=abfload('d:\data01.abf','info') 
% abfload will not load any data but return detailed information (header
% parameters) on the file in output variable h. d and si will be empty.
% In all other cases abfload will load data. Optional input parameters
% listed below (= all except the file name) must be specified as
% parameter/value pairs, e.g. as in 
%         d=abfload('d:\data01.abf','start',100,'stop','e');
%
% >>> INPUT VARIABLES >>>
% NAME        TYPE, DEFAULT      DESCRIPTION
% fn          char array         abf data file name
% start       scalar, 0          only gap-free-data: start of cutout to be 
%                                 read (unit: s)
% stop        scalar or char,    only gap-free-data: end of cutout to be  
%             'e'                 read (unit: sec). May be set to 'e' (end 
%                                 of file).
% sweeps      1d-array or char,  only episodic data: sweep numbers to be 
%             'a'                 read. By default, all sweeps will be read
%                                 ('a').
% channels    cell array         names of channels to be read, like 
%              or char, 'a'       {'IN 0','IN 8'} (make sure spelling is
%                                 100% correct, including blanks). If set 
%                                 to 'a', all channels will be read. 
%                                 *****************************************
%                                 NOTE: channel order in output variable d
%                                 ignores the order in 'channels', and
%                                 instead always matches the order inherent
%                                 to the abf file, to be retrieved in
%                                 output variable h!
%                                 *****************************************
% chunk       scalar, 0.05       only gap-free-data: the elementary chunk  
%                                 size (megabytes) to be used for the 
%                                 'discontinuous' mode of reading data 
%                                 (fewer channels to be read than exist)
% machineF    char array,        the 'machineformat' input parameter of the
%              'ieee-le'          matlab fopen function. 'ieee-le' is the 
%                                 correct option for windows; depending on 
%                                 the platform the data were recorded/shall
%                                 be read by abfload 'ieee-be' is the 
%                                 alternative.
% << OUTPUT VARIABLES <<<
% NAME  TYPE            DESCRIPTION
% data                     the data read, the format depending on the record-
%                        ing mode
%   1. GAP-FREE:
%   2d array        2d array of size 
%                    <data pts> by <number of chans>
%                    Examples of access:
%                    d(:,2)       data from channel 2 at full length
%                    d(1:100,:)   first 100 data points from all channels
%   2. EPISODIC FIXED-LENGTH/WAVEFORM FIXED-LENGTH/HIGH-SPEED OSCILLOSCOPE:
%   3d array        3d array of size 
%                    <data pts per sweep> by <number of chans> by <number 
%                    of sweeps>.
%                    Examples of access:
%                    d(:,2,:)            a matrix containing all episodes 
%                                        (at full length) of the second 
%                                        channel in its columns
%                    d(1:200,:,[1 11])   contains first 200 data points of 
%                                        episodes 1 and 11 of all channels
%   3. EPISODIC VARIABLE-LENGTH:
%   cell array      cell array whose elements correspond to single sweeps. 
%                    Each element is a (regular) array of size
%                    <data pts per sweep> by <number of chans>
%                    Examples of access:
%                    d{1}            a 2d-array which contains episode 1 
%                                    (all of it, all channels)
%                    d{2}(1:100,2)   a 1d-array containing the first 100
%                                    data points of channel 2 in episode 1
% si    scalar           the sampling interval in us
% fhead     struct           information on file (selected header parameters)
% 
% CONTRIBUTORS
%   Original version by Harald Hentschke (harald.hentschke@uni-tuebingen.de)
%   Extended to abf version 2.0 by Forrest Collman (fcollman@Princeton.edu)
%   pvpmod.m by Ulrich Egert (egert@bccn.uni-freiburg.de)

% -------------------------------------------------------------------------
%                       PART 1: check of input vars
% -------------------------------------------------------------------------
disp(['** ' mfilename])
% --- defaults   
% gap-free
start=0.0;
stop='e';
% episodic
sweeps='a';
% general
channels='a';
% the size of data chunks (see above) in Mb. 0.05 Mb is an empirical value
% which works well for abf with 6-16 channels and recording durations of 
% 5-30 min
chunk=0.05;
machineF='ieee-le';
verbose=1;
% if first and only optional input argument is string 'info' the user's
% request is to obtain information on the file (header parameters), so set
% flag accordingly
if nargin==2 && ischar(varargin{1}) && strcmp('info',varargin{1})
  doLoadData=false;
else
  doLoadData=true;
  % assign values of optional input parameters if any were given
  pvpmod(varargin);
end

% some constants
BLOCKSIZE=512;
% output variables
data=[]; 
si=[];
fhead=[];
if ischar(stop)
  if ~strcmpi(stop,'e')
    error('input parameter ''stop'' must be specified as ''e'' (=end of recording) or as a scalar');
  end
end
% check existence of file
if ~exist(fn,'file')
  error(['could not find file ' fn]);
end

% -------------------------------------------------------------------------
%                       PART 2a: determine abf version
% -------------------------------------------------------------------------
disp(['opening ' fn '...']);
[fid,messg]=fopen(fn,'r',machineF);
if fid == -1
  error(messg);
end
% on the occasion, determine absolute file size
fseek(fid,0,'eof');
fileSz=ftell(fid);
fseek(fid,0,'bof');

% *** read value of parameter 'fFileSignature' (i.e. abf version) from header ***
sz=4;
[fFileSignature,n]=fread(fid,sz,'uchar=>char');
if n~=sz
  fclose(fid);
  error('something went wrong reading value(s) for fFileSignature');
end
% rewind
fseek(fid,0,'bof');
% transpose
fFileSignature=fFileSignature';

% one of the first checks must be whether file signature is valid
switch fFileSignature
  case 'ABF ' % ** note the blank
    % ************************
    %     abf version < 2.0
    % ************************
  case 'ABF2'
    % ************************
    %     abf version >= 2.0
    % ************************
  otherwise
    error('unknown or incompatible file signature');
end

% -------------------------------------------------------------------------
%    PART 2b: define file information ('header' parameters) of interest
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


% define header proper depending on ABF version by call to local function 
headPar=define_header(fFileSignature);
% define all sections that there are
Sections=define_Sections;
% define a few of these (currently, only the TagInfo section is used for
% all versions of ABF, but that may change in the future)
ProtocolInfo=define_ProtocolInfo;
ADCInfo=define_ADCInfo;
TagInfo=define_TagInfo;
UserListInfo=define_UserListInfo;
EpochPerDACInfo=define_EpochPerDACInfo;

% -------------------------------------------------------------------------
%    PART 2c: read parameters of interest
% -------------------------------------------------------------------------
% convert headPar to struct
s=cell2struct(headPar,{'name','offs','numType','value'},2);
numOfParams=size(s,1);
clear tmp headPar;

% convert names in structure to variables and read value from header
for g=1:numOfParams
  if fseek(fid, s(g).offs,'bof')~=0
    fclose(fid);
    error(['something went wrong locating ' s(g).name]);
  end
  sz=length(s(g).value);
  % use dynamic field names
  [fhead.(s(g).name),n]=fread(fid,sz,s(g).numType);
  if n~=sz
    fclose(fid);
    error(['something went wrong reading value(s) for ' s(g).name]);
  end
end
% file signature needs to be transposed
fhead.fFileSignature=fhead.fFileSignature';
% several header parameters need a fix or version-specific refinement:
if strcmp(fhead.fFileSignature,'ABF2')
  % h.fFileVersionNumber needs to be converted from an array of integers to
  % a float
  fhead.fFileVersionNumber=fhead.fFileVersionNumber(4)+fhead.fFileVersionNumber(3)*.1...
    +fhead.fFileVersionNumber(2)*.001+fhead.fFileVersionNumber(1)*.0001;
  % convert ms to s
  fhead.lFileStartTime=fhead.uFileStartTimeMS*.001;
else
  % h.fFileVersionNumber is a float32 the value of which is sometimes a
  % little less than what it should be (e.g. 1.6499999 instead of 1.65)
  fhead.fFileVersionNumber=.001*round(fhead.fFileVersionNumber*1000);
  % in abf < 2.0 two parameters are needed to obtain the file start time
  % with millisecond precision - let's integrate both into parameter
  % lFileStartTime (unit: s) so that nFileStartMillisecs will not be needed
  fhead.lFileStartTime=fhead.lFileStartTime+fhead.nFileStartMillisecs*.001;
end

if fhead.fFileVersionNumber>=2
  % -----------------------------------------------------------------------
  % *** read file information that has moved from the header section to
  % other sections in ABF version >= 2.0 and assign selected values to
  % fields of 'generic' header variable h ***
  % -----------------------------------------------------------------------
  % --- read in the Sections
  Sects=cell2struct(Sections,{'name'},2);
  numOfSections=length(Sections);
  offset=76;
  % this creates all sections (ADCSection, ProtocolSection, etc.)
  for i=1:numOfSections
    eval([Sects(i).name '=ReadSectionInfo(fid,offset);']);
    offset=offset+4+4+8;
  end
  % --- read in the StringsSection and use some fields (to retrieve
  % information on the names of recorded channels and the units as well as
  % the recording protocol used)
  fseek(fid,StringsSection.uBlockIndex*BLOCKSIZE,'bof');
  BigString=fread(fid,StringsSection.uBytes,'char');
  % extract path and name of protocol file (excluding drive letter):
  % - find the first occurrence of a backslash in BigString
  tmpIx1=strfind(lower(char(BigString)'),'\');
  % - find '.pro' as this is the file extension of protocol files
  tmpIx2=strfind(lower(char(BigString)'),'.pro');
  % - extract everything in between and place in field of header struct h
  if ~isempty(tmpIx1) && ~isempty(tmpIx2) 
    fhead.protocolName=char(BigString(tmpIx1:tmpIx2(1)+3))';
  else
    fhead.protocolName='protocol name could not be identified';
  end
  % this is a hack: determine where either of strings 'clampex',
  % 'clampfit', 'axoscope' or patchxpress' begin
  progString={'clampex','clampfit','axoscope','patchxpress'};
  goodstart=[];
  for i=1:numel(progString)
    goodstart=cat(1,goodstart,strfind(lower(char(BigString)'),progString{i}));
  end
  if isempty(goodstart)
    error('problems in StringsSection: unrecognized acquisition program');
  end
  BigString=BigString(goodstart(1):end)';
  stringends=find(BigString==0);
  stringends=[0 stringends];
  Strings=cell(1,length(stringends)-1);
  for i=1:length(stringends)-1
    Strings{i}=char(BigString(stringends(i)+1:stringends(i+1)-1));
  end
  fhead.recChNames=[];
  fhead.recChUnits=[];
  
  % --- read in the ADCSection & copy some values to header h
  for i=1:ADCSection.llNumEntries
    ADCsec(i)=ReadSection(fid,ADCSection.uBlockIndex*BLOCKSIZE+ADCSection.uBytes*(i-1),ADCInfo);
    ii=ADCsec(i).nADCNum+1;
    fhead.nADCSamplingSeq(i)=ADCsec(i).nADCNum;
    fhead.recChNames=char(fhead.recChNames, Strings{ADCsec(i).lADCChannelNameIndex});
    unitsIndex=ADCsec(i).lADCUnitsIndex;
    if unitsIndex>0
        fhead.recChUnits=char(fhead.recChUnits, Strings{ADCsec(i).lADCUnitsIndex});
    else
        fhead.recChUnits=char(fhead.recChUnits,'nil');
    end
    fhead.nTelegraphEnable(ii)=ADCsec(i).nTelegraphEnable;
    fhead.fTelegraphAdditGain(ii)=ADCsec(i).fTelegraphAdditGain;
    fhead.fInstrumentScaleFactor(ii)=ADCsec(i).fInstrumentScaleFactor;
    fhead.fSignalGain(ii)=ADCsec(i).fSignalGain;
    fhead.fADCProgrammableGain(ii)=ADCsec(i).fADCProgrammableGain;
    fhead.fInstrumentOffset(ii)=ADCsec(i).fInstrumentOffset;
    fhead.fSignalOffset(ii)=ADCsec(i).fSignalOffset;
  end
  % --- read in the protocol section & copy some values to header h
  ProtocolSec=ReadSection(fid,ProtocolSection.uBlockIndex*BLOCKSIZE,ProtocolInfo);
  fhead.nOperationMode=ProtocolSec.nOperationMode;
  fhead.fSynchTimeUnit=ProtocolSec.fSynchTimeUnit;
  
  fhead.nADCNumChannels=ADCSection.llNumEntries;
  fhead.lActualAcqLength=DataSection.llNumEntries;
  fhead.lDataSectionPtr=DataSection.uBlockIndex;
  fhead.nNumPointsIgnored=0;
  % in ABF version < 2.0 h.fADCSampleInterval is the sampling interval
  % defined as
  %     1/(sampling freq*number_of_channels)
  % so divide ProtocolSec.fADCSequenceInterval by the number of
  % channels
  fhead.fADCSampleInterval=ProtocolSec.fADCSequenceInterval/fhead.nADCNumChannels;
  fhead.fADCRange=ProtocolSec.fADCRange;
  fhead.lADCResolution=ProtocolSec.lADCResolution;
  
  % --- read in the Epoch section and copy some values to header h
  for i=1:EpochPerDACSection.llNumEntries
    EPDsec(i)=ReadSection(fid,EpochPerDACSection.uBlockIndex*BLOCKSIZE+EpochPerDACSection.uBytes*(i-1),EpochPerDACInfo);
    ii=EPDsec(i).nEpochNum+1;
    fhead.nEpochNum(ii)=EPDsec(i).nEpochNum;
    fhead.nDACNum(ii)=EPDsec(i).nDACNum;
    fhead.nEpochType(ii)=EPDsec(i).nEpochType;
    fhead.fEpochInitLevel(ii)=EPDsec(i).fEpochInitLevel;
    fhead.fEpochLevelInc(ii)=EPDsec(i).fEpochLevelInc;
    fhead.lEpochInitDuration(ii)=EPDsec(i).lEpochInitDuration;
    fhead.lEpochDurationInc(ii)=EPDsec(i).lEpochDurationInc;
    fhead.lEpochPulsePeriod(ii)=EPDsec(i).lEpochPulsePeriod;
    fhead.lEpochPulseWidth(ii)=EPDsec(i).lEpochPulseWidth;
  end
  % --- read in the user list section - unclear how to get the values
  % (commented code lines below are a first guess)
  UserListSec=ReadSection(fid,UserListSection.uBlockIndex*BLOCKSIZE,UserListInfo);
  %   fseek(fid,UserListSec.lULParamValueListIndex*BLOCKSIZE,'bof');
  %   nada=fread(fid,100,'uchar=>char')
  
  % --- in contrast to procedures with all other sections do not read the 
  % sync array section but rather copy the values of its fields to the
  % corresponding fields of h
  fhead.lSynchArrayPtr=SynchArraySection.uBlockIndex;
  fhead.lSynchArraySize=SynchArraySection.llNumEntries;
else
  % -------------------------------------------------------------------------
  % *** here, do the inverse: in ABF version<2 files extract information
  % from header variable h and place it in corresponding new section
  % variable(s)
  % -------------------------------------------------------------------------
  TagSection.llNumEntries=fhead.lNumTagEntries;
  TagSection.uBlockIndex=fhead.lTagSectionPtr;
  TagSection.uBytes=64;
end
% -------------------------------------------------------------------------
%    PART 2d: groom parameters & perform some plausibility checks
% -------------------------------------------------------------------------
if fhead.lActualAcqLength<fhead.nADCNumChannels
  fclose(fid);
  error('less data points than sampled channels in file');
end
% the numerical value of all recorded channels (numbers 0..15)
recChIdx=fhead.nADCSamplingSeq(1:fhead.nADCNumChannels);
% the corresponding indices into loaded data d
recChInd=1:length(recChIdx);
if fhead.fFileVersionNumber<2
  % the channel names, e.g. 'IN 8' (for ABF version 2.0 these have been
  % extracted above at this point)
  fhead.recChNames=(reshape(char(fhead.sADCChannelName),10,16))';
  fhead.recChNames=fhead.recChNames(recChIdx+1,:);
  % same with signal units
  fhead.recChUnits=(reshape(char(fhead.sADCUnits),8,16))';
  fhead.recChUnits=fhead.recChUnits(recChIdx+1,:);
end
% convert to cell arrays
fhead.recChNames=deblank(cellstr(fhead.recChNames));
fhead.recChUnits=deblank(cellstr(fhead.recChUnits));

% check whether requested channels exist
chInd=[];
eflag=0;
if ischar(channels)
  if strcmp(channels,'a')
    chInd=recChInd;
  else
    fclose(fid);
    error('input parameter ''channels'' must either be a cell array holding channel names or the single character ''a'' (=all channels)');
  end
else
  % check for requested channels which do not exist
  missingChan=setdiff(channels,fhead.recChNames);
  % identify requested channels among available ones
  [nil,chInd]=intersect(fhead.recChNames,channels);
  % ** index chInd must be sorted because intersect sorts h.recChNames
  % alphanumerically, which needs not necessarily correspond to the order
  % inherent in the abf file (e.g. if channels are named 'Lynx1 ... Lynx10
  % etc.)
  chInd=sort(chInd);
  if isempty(chInd) || ~isempty(missingChan)
    % set error flag to 1
    eflag=1;
  end
end
if eflag
  fclose(fid);
  disp('**** available channels:');
  disp(fhead.recChNames);
  disp(' ');
  disp('**** requested channels:');
  disp(channels);
  error('at least one of the requested channels does not exist in data file (see above)');
end
% display available channels if in info mode
if ~doLoadData
  disp('**** available channels:');
  disp(fhead.recChNames);
end

% gain of telegraphed instruments, if any
if fhead.fFileVersionNumber>=1.65
  addGain=fhead.nTelegraphEnable.*fhead.fTelegraphAdditGain;
  addGain(addGain==0)=1;
else
  addGain=ones(size(fhead.fTelegraphAdditGain));
end

% determine offset at which data start
switch fhead.nDataFormat
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
headOffset=fhead.lDataSectionPtr*BLOCKSIZE+fhead.nNumPointsIgnored*dataSz;
% h.fADCSampleInterval is the TOTAL sampling interval
fhead.si=fhead.fADCSampleInterval*fhead.nADCNumChannels;
% assign same value to si, which is an output variable
si=fhead.si;
if ischar(sweeps) && sweeps=='a'
  nSweeps=fhead.lActualEpisodes;
  sweeps=1:fhead.lActualEpisodes;
else
  nSweeps=length(sweeps);
end

% determine time unit in synch array section
switch fhead.fSynchTimeUnit
  case 0  
    % time information in synch array section is in terms of ticks
    fhead.synchArrTimeBase=1;
  otherwise
    % time information in synch array section is in terms of usec
    fhead.synchArrTimeBase=fhead.fSynchTimeUnit;
end

% read in the TagSection, do a few computations & write to h.tags
fhead.tags=[];
for i=1:TagSection.llNumEntries
  tmp=ReadSection(fid,TagSection.uBlockIndex*BLOCKSIZE+TagSection.uBytes*(i-1),TagInfo);
  % time of tag entry from start of experiment in s (corresponding expisode
  % number, if applicable, will be determined later)
  fhead.tags(i).timeSinceRecStart=tmp.lTagTime*fhead.synchArrTimeBase/1e6;
  fhead.tags(i).comment=char(tmp.sComment)';
end

% -------------------------------------------------------------------------
%    PART 3: read data (note: from here on code is generic and abf version
%    should not matter)
% -------------------------------------------------------------------------
switch fhead.nOperationMode
  case 1
    disp('data were acquired in event-driven variable-length mode');
    if fhead.fFileVersionNumber>=2.0
      errordlg('abfload currently does not work with data acquired in event-driven variable-length mode and ABF version 2.0','ABF version issue');
    else
      if (fhead.lSynchArrayPtr<=0 || fhead.lSynchArraySize<=0)
        fclose(fid);
        error('internal variables ''lSynchArraynnn'' are zero or negative');
      end
      % the byte offset at which the SynchArraySection starts
      fhead.lSynchArrayPtrByte=BLOCKSIZE*fhead.lSynchArrayPtr;
      % before reading Synch Arr parameters check if file is big enough to hold them
      % 4 bytes/long, 2 values per episode (start and length)
      if fhead.lSynchArrayPtrByte+2*4*fhead.lSynchArraySize<fileSz
        fclose(fid);
        error('file seems not to contain complete Synch Array Section');
      end
      if fseek(fid,fhead.lSynchArrayPtrByte,'bof')~=0
        fclose(fid);
        error('something went wrong positioning file pointer to Synch Array Section');
      end
      [synchArr,n]=fread(fid,fhead.lSynchArraySize*2,'int32');
      if n~=fhead.lSynchArraySize*2
        fclose(fid);
        error('something went wrong reading synch array section');
      end
      % make synchArr a h.lSynchArraySize x 2 matrix
      synchArr=permute(reshape(synchArr',2,fhead.lSynchArraySize),[2 1]);
      % the length of episodes in sample points
      segLengthInPts=synchArr(:,2)/fhead.synchArrTimeBase;
      % the starting ticks of episodes in sample points WITHIN THE DATA FILE
      segStartInPts=cumsum([0 (segLengthInPts(1:end-1))']*dataSz)+headOffset;
      % start time (synchArr(:,1)) has to be divided by h.nADCNumChannels to get true value
      % go to data portion
      if fseek(fid,headOffset,'bof')~=0
        fclose(fid);
        error('something went wrong positioning file pointer (too few data points ?)');
      end
      % ** load data if requested
      if doLoadData
        data=cell(1,nSweeps);
        for i=1:nSweeps
          % if selected sweeps are to be read, seek correct position
          if ~isequal(nSweeps,fhead.lActualEpisodes)
            fseek(fid,segStartInPts(sweeps(i)),'bof');
          end
          [tmpd,n]=fread(fid,segLengthInPts(sweeps(i)),precision);
          if n~=segLengthInPts(sweeps(i))
            warning(['something went wrong reading episode ' int2str(sweeps(i)) ': ' segLengthInPts(sweeps(i)) ' points should have been read, ' int2str(n) ' points actually read']);
          end
          fhead.dataPtsPerChan=n/fhead.nADCNumChannels;
          if rem(n,fhead.nADCNumChannels)>0
            fclose(fid);
            error('number of data points in episode not OK');
          end
          % separate channels..
          tmpd=reshape(tmpd,fhead.nADCNumChannels,fhead.dataPtsPerChan);
          % retain only requested channels
          tmpd=tmpd(chInd,:);
          tmpd=tmpd';
          % if data format is integer, scale appropriately; if it's float, tmpd is fine
          if ~fhead.nDataFormat
            for j=1:length(chInd)
              ch=recChIdx(chInd(j))+1;
              tmpd(:,j)=tmpd(:,j)/(fhead.fInstrumentScaleFactor(ch)*fhead.fSignalGain(ch)*fhead.fADCProgrammableGain(ch)*addGain(ch))...
                *fhead.fADCRange/fhead.lADCResolution+fhead.fInstrumentOffset(ch)-fhead.fSignalOffset(ch);
            end
          end
          % now place in cell array, an element consisting of one sweep with channels in columns
          data{i}=tmpd;
        end
      end
    end
    
  case {2,4,5}
    if fhead.nOperationMode==2
      disp('data were acquired in event-driven fixed-length mode');
    elseif fhead.nOperationMode==4
      disp('data were acquired in high-speed oscilloscope mode');
    else
      disp('data were acquired in waveform fixed-length mode');
    end
    % extract timing information on sweeps
    if (fhead.lSynchArrayPtr<=0 || fhead.lSynchArraySize<=0)
      fclose(fid);
      error('internal variables ''lSynchArraynnn'' are zero or negative');
    end
    % the byte offset at which the SynchArraySection starts
    fhead.lSynchArrayPtrByte=BLOCKSIZE*fhead.lSynchArrayPtr;
    % before reading Synch Arr parameters check if file is big enough to hold them
    % 4 bytes/long, 2 values per episode (start and length)
    if fhead.lSynchArrayPtrByte+2*4*fhead.lSynchArraySize>fileSz
      fclose(fid);
      error('file seems not to contain complete Synch Array Section');
    end
    if fseek(fid,fhead.lSynchArrayPtrByte,'bof')~=0
      fclose(fid);
      error('something went wrong positioning file pointer to Synch Array Section');
    end
    [synchArr,n]=fread(fid,fhead.lSynchArraySize*2,'int32');
    if n~=fhead.lSynchArraySize*2
      fclose(fid);
      error('something went wrong reading synch array section');
    end
    % make synchArr a h.lSynchArraySize x 2 matrix
    synchArr=permute(reshape(synchArr',2,fhead.lSynchArraySize),[2 1]);
    if numel(unique(synchArr(:,2)))>1
      fclose(fid);
      error('sweeps of unequal length in file recorded in fixed-length mode');
    end
    % the length of sweeps in sample points (**note: parameter lLength of
    % the ABF synch section is expressed in samples (ticks) whereas
    % parameter lStart is given in synchArrTimeBase units)
    fhead.sweepLengthInPts=synchArr(1,2)/fhead.nADCNumChannels;
    % the starting ticks of episodes in sample points (t0=1=beginning of
    % recording)
    fhead.sweepStartInPts=synchArr(:,1)*(fhead.synchArrTimeBase/fhead.fADCSampleInterval/fhead.nADCNumChannels);
    % recording start and stop times in seconds from midnight
    fhead.recTime=fhead.lFileStartTime;
    fhead.recTime=fhead.recTime+[0  (1e-6*(fhead.sweepStartInPts(end)+fhead.sweepLengthInPts))*fhead.fADCSampleInterval*fhead.nADCNumChannels];
    % determine first point and number of points to be read
    startPt=0;
    fhead.dataPts=fhead.lActualAcqLength;
    fhead.dataPtsPerChan=fhead.dataPts/fhead.nADCNumChannels;
    if rem(fhead.dataPts,fhead.nADCNumChannels)>0 || rem(fhead.dataPtsPerChan,fhead.lActualEpisodes)>0
      fclose(fid);
      error('number of data points not OK');
    end
    % temporary helper var
    dataPtsPerSweep=fhead.sweepLengthInPts*fhead.nADCNumChannels;
    if fseek(fid,startPt*dataSz+headOffset,'bof')~=0
      fclose(fid);
      error('something went wrong positioning file pointer (too few data points ?)');
    end
    % the starting ticks of episodes in sample points WITHIN THE DATA FILE
    selectedSegStartInPts=((sweeps-1)*dataPtsPerSweep)*dataSz+headOffset;
    % ** load data if requested
    if doLoadData
      % preallocate d
      data=zeros(fhead.sweepLengthInPts,length(chInd),nSweeps);
      for i=1:nSweeps
        status=fseek(fid,selectedSegStartInPts(i),'bof');
        if status==-1
          fclose(fid);
          error(['something went wrong reading episode ' int2str(sweeps(i)) '; file pointer beyond file limits (check sweeps)']);
        end
        [tmpd,n]=fread(fid,dataPtsPerSweep,precision);
        if n~=dataPtsPerSweep
          fclose(fid);
          error(['something went wrong reading episode ' int2str(sweeps(i)) ': ' dataPtsPerSweep ' points should have been read, ' int2str(n) ' points actually read']);
        end
        fhead.dataPtsPerChan=n/fhead.nADCNumChannels;
        if rem(n,fhead.nADCNumChannels)>0
          fclose(fid);
          error('number of data points in episode not OK');
        end
        % separate channels..
        tmpd=reshape(tmpd,fhead.nADCNumChannels,fhead.dataPtsPerChan);
        % retain only requested channels
        tmpd=tmpd(chInd,:);
        tmpd=tmpd';
        % if data format is integer, scale appropriately; if it's float, d is fine
        if ~fhead.nDataFormat
          for j=1:length(chInd)
            ch=recChIdx(chInd(j))+1;
            tmpd(:,j)=tmpd(:,j)/(fhead.fInstrumentScaleFactor(ch)*fhead.fSignalGain(ch)*fhead.fADCProgrammableGain(ch)*addGain(ch))...
              *fhead.fADCRange/fhead.lADCResolution+fhead.fInstrumentOffset(ch)-fhead.fSignalOffset(ch);
          end
        end
        % now fill 3d array
        data(:,:,i)=tmpd;
      end
    end
    
  case 3
    disp('data were acquired in gap-free mode');
    % from start, stop, headOffset and h.fADCSampleInterval calculate first point to be read
    %  and - unless stop is given as 'e' - number of points
    startPt=floor(1e6*start*(1/fhead.fADCSampleInterval));
    % this corrects undesired shifts in the reading frame due to rounding errors in the previous calculation
    startPt=floor(startPt/fhead.nADCNumChannels)*fhead.nADCNumChannels;
    % if stop is a char array, it can only be 'e' at this point (other values would have
    % been caught above)
    if ischar(stop)
      fhead.dataPtsPerChan=fhead.lActualAcqLength/fhead.nADCNumChannels-floor(1e6*start/fhead.si);
      fhead.dataPts=fhead.dataPtsPerChan*fhead.nADCNumChannels;
    else
      fhead.dataPtsPerChan=floor(1e6*(stop-start)*(1/fhead.si));
      fhead.dataPts=fhead.dataPtsPerChan*fhead.nADCNumChannels;
      if fhead.dataPts<=0
        fclose(fid);
        error('start is larger than or equal to stop');
      end
    end
    if rem(fhead.dataPts,fhead.nADCNumChannels)>0
      fclose(fid);
      error('number of data points not OK');
    end
    tmp=1e-6*fhead.lActualAcqLength*fhead.fADCSampleInterval;
    if verbose
      disp(['total length of recording: ' num2str(tmp,'%5.1f') ' s ~ ' num2str(tmp/60,'%3.0f') ' min']);
      disp(['sampling interval: ' num2str(fhead.si,'%5.0f') ' µs']);
      % 8 bytes per data point expressed in Mb
      disp(['memory requirement for complete upload in matlab: '...
        num2str(round(8*fhead.lActualAcqLength/2^20)) ' MB']);
    end
    % recording start and stop times in seconds from midnight
    fhead.recTime=fhead.lFileStartTime;
    fhead.recTime=[fhead.recTime fhead.recTime+tmp];
    if fseek(fid,startPt*dataSz+headOffset,'bof')~=0
      fclose(fid);
      error('something went wrong positioning file pointer (too few data points ?)');
    end
    if doLoadData
      % *** decide on the most efficient way to read data:
      % (i) all (of one or several) channels requested: read, done
      % (ii) one (of several) channels requested: use the 'skip' feature of
      % fread
      % (iii) more than one but not all (of several) channels requested:
      % 'discontinuous' mode of reading data. Read a reasonable chunk of data
      % (all channels), separate channels, discard non-requested ones (if
      % any), place data in preallocated array, repeat until done. This is
      % faster than reading the data in one big lump, separating channels and
      % discarding the ones not requested
      if length(chInd)==1 && fhead.nADCNumChannels>1
        % --- situation (ii)
        % jump to proper reading frame position in file
        if fseek(fid,(chInd-1)*dataSz,'cof')~=0
          fclose(fid);
          error('something went wrong positioning file pointer (too few data points ?)');
        end
        % read, skipping h.nADCNumChannels-1 data points after each read
        [data,n]=fread(fid,fhead.dataPtsPerChan,precision,dataSz*(fhead.nADCNumChannels-1));
        if n~=fhead.dataPtsPerChan
          fclose(fid);
          error(['something went wrong reading file (' int2str(fhead.dataPtsPerChan) ' points should have been read, ' int2str(n) ' points actually read']);
        end
      elseif length(chInd)/fhead.nADCNumChannels<1
        % --- situation (iii)
        % prepare chunkwise upload:
        % preallocate d
        data=nan(fhead.dataPtsPerChan,length(chInd));
        % the number of data points corresponding to the maximal chunk size,
        % rounded off such that from each channel the same number of points is
        % read (do not forget that each data point will by default be made a
        % double of 8 bytes, no matter what the original data format is)
        chunkPtsPerChan=floor(chunk*2^20/8/fhead.nADCNumChannels);
        chunkPts=chunkPtsPerChan*fhead.nADCNumChannels;
        % the number of those chunks..
        nChunk=floor(fhead.dataPts/chunkPts);
        % ..and the remainder
        restPts=fhead.dataPts-nChunk*chunkPts;
        restPtsPerChan=restPts/fhead.nADCNumChannels;
        % chunkwise row indices into d
        dix=(1:chunkPtsPerChan:fhead.dataPtsPerChan)';
        dix(:,2)=dix(:,1)+chunkPtsPerChan-1;
        dix(end,2)=fhead.dataPtsPerChan;
        if verbose && nChunk
          disp(['reading file in ' int2str(nChunk) ' chunks of ~' num2str(chunk) ' Mb']);
        end
        % do it: if no remainder exists loop through all rows of dix,
        % otherwise spare last row for the lines below (starting with
        % 'if restPts')
        for ci=1:size(dix,1)-(restPts>0)
          [tmpd,n]=fread(fid,chunkPts,precision);
          if n~=chunkPts
            fclose(fid);
            error(['something went wrong reading chunk #' int2str(ci) ' (' ...
              int2str(chunkPts) ' points should have been read, ' int2str(n) ' points actually read']);
          end
          % separate channels..
          tmpd=reshape(tmpd,fhead.nADCNumChannels,chunkPtsPerChan);
          data(dix(ci,1):dix(ci,2),:)=tmpd(chInd,:)';
        end
        % collect the rest, if any
        if restPts
          [tmpd,n]=fread(fid,restPts,precision);
          if n~=restPts
            fclose(fid);
            error(['something went wrong reading last chunk (' ...
              int2str(restPts) ' points should have been read, ' int2str(n) ' points actually read']);
          end
          % separate channels..
          tmpd=reshape(tmpd,fhead.nADCNumChannels,restPtsPerChan);
          data(dix(end,1):dix(end,2),:)=tmpd(chInd,:)';
        end
      else
        % --- situation (i)
        [data,n]=fread(fid,fhead.dataPts,precision);
        if n~=fhead.dataPts
          fclose(fid);
          error(['something went wrong reading file (' int2str(fhead.dataPts) ' points should have been read, ' int2str(n) ' points actually read']);
        end
        % separate channels..
        data=reshape(data,fhead.nADCNumChannels,fhead.dataPtsPerChan);
        data=data';
      end
      % if data format is integer, scale appropriately; if it's float, d is fine
      if ~fhead.nDataFormat
        for j=1:length(chInd)
          ch=recChIdx(chInd(j))+1;
          data(:,j)=data(:,j)/(fhead.fInstrumentScaleFactor(ch)*fhead.fSignalGain(ch)*fhead.fADCProgrammableGain(ch)*addGain(ch))...
            *fhead.fADCRange/fhead.lADCResolution+fhead.fInstrumentOffset(ch)-fhead.fSignalOffset(ch);
        end
      end
    end
  otherwise
    disp('unknown recording mode -- returning empty matrix');
    data=[];
    fhead.si=[];
end
fclose(fid);

% finally, possibly add information on episode number to tags
if ~isempty(fhead.tags) && isfield(fhead,'sweepStartInPts')
  for i=1:numel(fhead.tags)
    tmp=find(fhead.tags(i).timeSinceRecStart>=fhead.sweepStartInPts/1e6*fhead.si);
    fhead.tags(i).episodeIndex=tmp(end);
  end
end


% ########################################################################
%                         LOCAL FUNCTIONS
% ########################################################################

function headPar=define_header(fileSig)
switch fileSig
 case 'ABF ' % ** note the blank
   % ************************
   %     abf version < 2.0
   % ************************
   %
   % temporary initializing var
   tmp=repmat(-1,1,16);
   % define vital header parameters and initialize them with -1: set up a
   % cell array (and convert it to a struct later on, which is more
   % convenient)
   % column order is
   %    name, position in header in bytes, type, value)
   headPar={
     'fFileSignature',0,'*char',[-1 -1 -1 -1];
     'fFileVersionNumber',4,'float32',-1;
     'nOperationMode',8,'int16',-1;
     'lActualAcqLength',10,'int32',-1;
     'nNumPointsIgnored',14,'int16',-1;
     'lActualEpisodes',16,'int32',-1;
     'lFileStartTime',24,'int32',-1;
     'lDataSectionPtr',40,'int32',-1;
     'lTagSectionPtr',44,'int32',-1;
     'lNumTagEntries',48,'int32',-1;
     'lSynchArrayPtr',92,'int32',-1;
     'lSynchArraySize',96,'int32',-1;
     'nDataFormat',100,'int16',-1;
     'nADCNumChannels', 120, 'int16', -1;
     'fADCSampleInterval',122,'float', -1;
     'fSynchTimeUnit',130,'float',-1;
     'lNumSamplesPerEpisode',138,'int32',-1;
     'lPreTriggerSamples',142,'int32',-1;
     'lEpisodesPerRun',146,'int32',-1;
     'fADCRange', 244, 'float', -1;
     'lADCResolution', 252, 'int32', -1;
     'nFileStartMillisecs', 366, 'int16', -1;
     'nADCPtoLChannelMap', 378, 'int16', tmp;
     'nADCSamplingSeq', 410, 'int16',  tmp;
     'sADCChannelName',442, 'uchar', repmat(tmp,1,10);
     'sADCUnits',602, 'uchar', repmat(tmp,1,8);
     'fADCProgrammableGain', 730, 'float', tmp;
     'fInstrumentScaleFactor', 922, 'float', tmp;
     'fInstrumentOffset', 986, 'float', tmp;
     'fSignalGain', 1050, 'float', tmp;
     'fSignalOffset', 1114, 'float', tmp;
     'nTelegraphEnable',4512,'int16',tmp;
     'fTelegraphAdditGain',4576,'float',tmp
     };
 case 'ABF2'
   % ************************
   %     abf version >= 2.0
   % ************************
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
end
function EpochPerDACInfo=define_EpochPerDACInfo 
EpochPerDACInfo={ 
    'nEpochNum','int16',1; 
    'nDACNum','int16',1; 
    'nEpochType','int16',1; 
    'fEpochInitLevel','float',1; 
    'fEpochLevelInc','float',1; 
    'lEpochInitDuration','int32',1; 
    'lEpochDurationInc','int32',1; 
    'lEpochPulsePeriod','int32',1; 
    'lEpochPulseWidth','int32',1; 
};
function Sections=define_Sections
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

function TagInfo=define_TagInfo
TagInfo={
   'lTagTime','int32',1;
   'sComment','char',56;
   'nTagType','int16',1;
   'nVoiceTagNumber_or_AnnotationIndex','int16',1;
};

function UserListInfo=define_UserListInfo
UserListInfo={
   'nListNum','int16',1;
   'nULEnable','int16',1;
   'nULParamToVary','int16',1;
   'nULRepeat','int16',1;
   'lULParamValueListIndex','int32',1;
   'sUnused','uchar',52;   
};

% // FUNCTION: ReadUserList
% // PURPOSE:  Reads the user list from the data file.
% //
% BOOL CABF2ProtocolReader::ReadUserList()
% {
%     MEMBERASSERT();
% 
%     BOOL bOK = TRUE;
%     if( m_FileInfo.UserListSection.uBlockIndex )
%     {
%         ABF_UserListInfo UserList;
%         ASSERT( m_FileInfo.UserListSection.uBytes == sizeof( UserList ) );
%         ASSERT( m_FileInfo.UserListSection.llNumEntries );
%         bOK &= m_pFI->Seek( LONGLONG(m_FileInfo.UserListSection.uBlockIndex) * ABF_BLOCKSIZE, FILE_BEGIN );
%         if( !bOK )
%             return FALSE;
% 
%         for( long i=0; i<m_FileInfo.UserListSection.llNumEntries; i++ )
%         {
%             bOK &= m_pFI->Read( &UserList, sizeof( UserList ) );
%             short u = UserList.nListNum;        
% 
%             m_pFH->nULEnable[u]      = 1;    
%             m_pFH->nULParamToVary[u] = UserList.nULParamToVary;          
%             m_pFH->nULRepeat[u]      = UserList.nULRepeat;
% 
%             bOK &= GetString( UserList.lULParamValueListIndex, m_pFH->sULParamValueList[u], ABF_USERLISTLEN );
%         }
%     }
%     return bOK;
% }

function Section=ReadSection(fid,offset,Format)
s=cell2struct(Format,{'name','numType','number'},2);
fseek(fid,offset,'bof');
for i=1:length(s)
 [Section.(s(i).name)]=fread(fid,s(i).number,s(i).numType);
end

function SectionInfo=ReadSectionInfo(fid,offset)
fseek(fid,offset,'bof');
SectionInfo.uBlockIndex=fread(fid,1,'uint32');
fseek(fid,offset+4,'bof');
SectionInfo.uBytes=fread(fid,1,'uint32');
fseek(fid,offset+8,'bof');
SectionInfo.llNumEntries=fread(fid,1,'int64');

function pvpmod(x)
% PVPMOD             - evaluate parameter/value pairs
% pvpmod(x) assigns the value x(i+1) to the parameter defined by the
% string x(i) in the calling workspace. This is useful to evaluate 
% <varargin> contents in an mfile, e.g. to change default settings 
% of any variable initialized before pvpmod(x) is called.
%
% (c) U. Egert 1998

% this loop is assigns the parameter/value pairs in x to the calling
% workspace.
if ~isempty(x)
  for i = 1:2:size(x,2)
     assignin('caller', x{i}, x{i+1});
  end
end



% 
% struct ABF_FileInfo
% {
%    UINT  uFileSignature;
%    UINT  uFileVersionNumber;
% 
%    // After this point there is no need to be the same as the ABF 1 equivalent.
%    UINT  uFileInfoSize;
% 
%    UINT  uActualEpisodes;
%    UINT  uFileStartDate;
%    UINT  uFileStartTimeMS;
%    UINT  uStopwatchTime;
%    short nFileType;
%    short nDataFormat;
%    short nSimultaneousScan;
%    short nCRCEnable;
%    UINT  uFileCRC;
%    GUID  FileGUID;
%    UINT  uCreatorVersion;
%    UINT  uCreatorNameIndex;
%    UINT  uModifierVersion;
%    UINT  uModifierNameIndex;
%    UINT  uProtocolPathIndex;   
% 
%    // New sections in ABF 2 - protocol stuff ...
%    ABF_Section ProtocolSection;           // the protocol
%    ABF_Section ADCSection;                // one for each ADC channel
%    ABF_Section DACSection;                // one for each DAC channel
%    ABF_Section EpochSection;              // one for each epoch
%    ABF_Section ADCPerDACSection;          // one for each ADC for each DAC
%    ABF_Section EpochPerDACSection;        // one for each epoch for each DAC
%    ABF_Section UserListSection;           // one for each user list
%    ABF_Section StatsRegionSection;        // one for each stats region
%    ABF_Section MathSection;
%    ABF_Section StringsSection;
% 
%    // ABF 1 sections ...
%    ABF_Section DataSection;            // Data
%    ABF_Section TagSection;             // Tags
%    ABF_Section ScopeSection;           // Scope config
%    ABF_Section DeltaSection;           // Deltas
%    ABF_Section VoiceTagSection;        // Voice Tags
%    ABF_Section SynchArraySection;      // Synch Array
%    ABF_Section AnnotationSection;      // Annotations
%    ABF_Section StatsSection;           // Stats config
%    
%    char  sUnused[148];     // size = 512 bytes
%    
%    ABF_FileInfo() 
%    { 
%       MEMSET_CTOR;
%       STATIC_ASSERT( sizeof( ABF_FileInfo ) == 512 );
% 
%       uFileSignature = ABF_FILESIGNATURE;
%       uFileInfoSize  = sizeof( ABF_FileInfo);
%    }
% 
% };
% 
% struct ABF_ProtocolInfo
% {
%    short nOperationMode;
%    float fADCSequenceInterval;
%    bool  bEnableFileCompression;
%    char  sUnused1[3];
%    UINT  uFileCompressionRatio;
% 
%    float fSynchTimeUnit;
%    float fSecondsPerRun;
%    long  lNumSamplesPerEpisode;
%    long  lPreTriggerSamples;
%    long  lEpisodesPerRun;
%    long  lRunsPerTrial;
%    long  lNumberOfTrials;
%    short nAveragingMode;
%    short nUndoRunCount;
%    short nFirstEpisodeInRun;
%    float fTriggerThreshold;
%    short nTriggerSource;
%    short nTriggerAction;
%    short nTriggerPolarity;
%    float fScopeOutputInterval;
%    float fEpisodeStartToStart;
%    float fRunStartToStart;
%    long  lAverageCount;
%    float fTrialStartToStart;
%    short nAutoTriggerStrategy;
%    float fFirstRunDelayS;
% 
%    short nChannelStatsStrategy;
%    long  lSamplesPerTrace;
%    long  lStartDisplayNum;
%    long  lFinishDisplayNum;
%    short nShowPNRawData;
%    float fStatisticsPeriod;
%    long  lStatisticsMeasurements;
%    short nStatisticsSaveStrategy;
% 
%    float fADCRange;
%    float fDACRange;
%    long  lADCResolution;
%    long  lDACResolution;
%    
%    short nExperimentType;
%    short nManualInfoStrategy;
%    short nCommentsEnable;
%    long  lFileCommentIndex;            
%    short nAutoAnalyseEnable;
%    short nSignalType;
% 
%    short nDigitalEnable;
%    short nActiveDACChannel;
%    short nDigitalHolding;
%    short nDigitalInterEpisode;
%    short nDigitalDACChannel;
%    short nDigitalTrainActiveLogic;
% 
%    short nStatsEnable;
%    short nStatisticsClearStrategy;
% 
%    short nLevelHysteresis;
%    long  lTimeHysteresis;
%    short nAllowExternalTags;
%    short nAverageAlgorithm;
%    float fAverageWeighting;
%    short nUndoPromptStrategy;
%    short nTrialTriggerSource;
%    short nStatisticsDisplayStrategy;
%    short nExternalTagType;
%    short nScopeTriggerOut;
% 
%    short nLTPType;
%    short nAlternateDACOutputState;
%    short nAlternateDigitalOutputState;
% 
%    float fCellID[3];
% 
%    short nDigitizerADCs;
%    short nDigitizerDACs;
%    short nDigitizerTotalDigitalOuts;
%    short nDigitizerSynchDigitalOuts;
%    short nDigitizerType;
% 
%    char  sUnused[304];     // size = 512 bytes
%    
%    ABF_ProtocolInfo() 
%    { 
%       MEMSET_CTOR; 
%       STATIC_ASSERT( sizeof( ABF_ProtocolInfo ) == 512 );
%    }
% };
% 
% struct ABF_MathInfo
% {
%    short nMathEnable;
%    short nMathExpression;
%    UINT  uMathOperatorIndex;     
%    UINT  uMathUnitsIndex;        
%    float fMathUpperLimit;
%    float fMathLowerLimit;
%    short nMathADCNum[2];
%    char  sUnused[16];
%    float fMathK[6];
% 
%    char  sUnused2[64];     // size = 128 bytes
%    
%    ABF_MathInfo()
%    { 
%       MEMSET_CTOR; 
%       STATIC_ASSERT( sizeof( ABF_MathInfo ) == 128 );
%    }
% };
% 
% struct ABF_ADCInfo
% {
%    // The ADC this struct is describing.
%    short nADCNum;
% 
%    short nTelegraphEnable;
%    short nTelegraphInstrument;
%    float fTelegraphAdditGain;
%    float fTelegraphFilter;
%    float fTelegraphMembraneCap;
%    short nTelegraphMode;
%    float fTelegraphAccessResistance;
% 
%    short nADCPtoLChannelMap;
%    short nADCSamplingSeq;
% 
%    float fADCProgrammableGain;
%    float fADCDisplayAmplification;
%    float fADCDisplayOffset;
%    float fInstrumentScaleFactor;
%    float fInstrumentOffset;
%    float fSignalGain;
%    float fSignalOffset;
%    float fSignalLowpassFilter;
%    float fSignalHighpassFilter;
% 
%    char  nLowpassFilterType;
%    char  nHighpassFilterType;
%    float fPostProcessLowpassFilter;
%    char  nPostProcessLowpassFilterType;
%    bool  bEnabledDuringPN;
% 
%    short nStatsChannelPolarity;
% 
%    long  lADCChannelNameIndex;
%    long  lADCUnitsIndex;
% 
%    char  sUnused[46];         // size = 128 bytes
%    
%    ABF_ADCInfo()
%    { 
%       MEMSET_CTOR; 
%       STATIC_ASSERT( sizeof( ABF_ADCInfo ) == 128 );
%    }
% };
% 
% struct ABF_DACInfo
% {
%    // The DAC this struct is describing.
%    short nDACNum;
% 
%    short nTelegraphDACScaleFactorEnable;
%    float fInstrumentHoldingLevel;
% 
%    float fDACScaleFactor;
%    float fDACHoldingLevel;
%    float fDACCalibrationFactor;
%    float fDACCalibrationOffset;
% 
%    long  lDACChannelNameIndex;
%    long  lDACChannelUnitsIndex;
% 
%    long  lDACFilePtr;
%    long  lDACFileNumEpisodes;
% 
%    short nWaveformEnable;
%    short nWaveformSource;
%    short nInterEpisodeLevel;
% 
%    float fDACFileScale;
%    float fDACFileOffset;
%    long  lDACFileEpisodeNum;
%    short nDACFileADCNum;
% 
%    short nConditEnable;
%    long  lConditNumPulses;
%    float fBaselineDuration;
%    float fBaselineLevel;
%    float fStepDuration;
%    float fStepLevel;
%    float fPostTrainPeriod;
%    float fPostTrainLevel;
%    short nMembTestEnable;
% 
%    short nLeakSubtractType;
%    short nPNPolarity;
%    float fPNHoldingLevel;
%    short nPNNumADCChannels;
%    short nPNPosition;
%    short nPNNumPulses;
%    float fPNSettlingTime;
%    float fPNInterpulse;
% 
%    short nLTPUsageOfDAC;
%    short nLTPPresynapticPulses;
% 
%    long  lDACFilePathIndex;
% 
%    float fMembTestPreSettlingTimeMS;
%    float fMembTestPostSettlingTimeMS;
% 
%    short nLeakSubtractADCIndex;
% 
%    char  sUnused[124];     // size = 256 bytes
%    
%    ABF_DACInfo()
%    { 
%       MEMSET_CTOR; 
%       STATIC_ASSERT( sizeof( ABF_DACInfo ) == 256 );
%    }
% };
% 
% struct ABF_EpochInfoPerDAC
% {
%    // The Epoch / DAC this struct is describing.
%    short nEpochNum;
%    short nDACNum;
% 
%    // One full set of epochs (ABF_EPOCHCOUNT) for each DAC channel ...
%    short nEpochType;
%    float fEpochInitLevel;
%    float fEpochLevelInc;
%    long  lEpochInitDuration;  
%    long  lEpochDurationInc;
%    long  lEpochPulsePeriod;
%    long  lEpochPulseWidth;
% 
%    char  sUnused[18];      // size = 48 bytes
%    
%    ABF_EpochInfoPerDAC()
%    { 
%       MEMSET_CTOR; 
%       STATIC_ASSERT( sizeof( ABF_EpochInfoPerDAC ) == 48 );
%    }
% };
% 
% struct ABF_EpochInfo
% {
%    // The Epoch this struct is describing.
%    short nEpochNum;
% 
%    // Describes one epoch
%    short nDigitalValue;
%    short nDigitalTrainValue;
%    short nAlternateDigitalValue;
%    short nAlternateDigitalTrainValue;
%    bool  bEpochCompression;   // Compress the data from this epoch using uFileCompressionRatio
% 
%    char  sUnused[21];      // size = 32 bytes
%    
%    ABF_EpochInfo()
%    { 
%       MEMSET_CTOR; 
%       STATIC_ASSERT( sizeof( ABF_EpochInfo ) == 32 );
%    }
% };
% 
% struct ABF_StatsRegionInfo
% { 
%    // The stats region this struct is describing.
%    short nRegionNum;
%    short nADCNum;
% 
%    short nStatsActiveChannels;
%    short nStatsSearchRegionFlags;
%    short nStatsSelectedRegion;
%    short nStatsSmoothing;
%    short nStatsSmoothingEnable;
%    short nStatsBaseline;
%    long  lStatsBaselineStart;
%    long  lStatsBaselineEnd;
% 
%    // Describes one stats region
%    long  lStatsMeasurements;
%    long  lStatsStart;
%    long  lStatsEnd;
%    short nRiseBottomPercentile;
%    short nRiseTopPercentile;
%    short nDecayBottomPercentile;
%    short nDecayTopPercentile;
%    short nStatsSearchMode;
%    short nStatsSearchDAC;
%    short nStatsBaselineDAC;
% 
%    char  sUnused[78];   // size = 128 bytes
%    
%    ABF_StatsRegionInfo()
%    { 
%       MEMSET_CTOR; 
%       STATIC_ASSERT( sizeof( ABF_StatsRegionInfo ) == 128 );
%    }
% };
% 
% struct ABF_UserListInfo
% {
%    // The user list this struct is describing.
%    short nListNum;
% 
%    // Describes one user list
%    short nULEnable;
%    short nULParamToVary;
%    short nULRepeat;
%    long  lULParamValueListIndex;
% 
%    char  sUnused[52];   // size = 64 bytes
%    
%    ABF_UserListInfo()
%    { 
%       MEMSET_CTOR; 
%       STATIC_ASSERT( sizeof( ABF_UserListInfo ) == 64 );
%    }
% };*/=