function Dirs = setgetDirs(varargin)
% If no Arguments are given, return all the globally defined directories, 
% plus the standard directory
% User supplied directories overwrite the one in GV and the Default Dirs.
%
% current Directories are: Speaker, STRF, TIT, TIIT, Results, Colormaps
%
% CALL: setgetDirs('STRF','C:\project\STRFs')
global Dirs;
if ~isempty(varargin) & iscell(varargin) & iscell(varargin{1}) varargin = varargin{1}; end
Sep = HF_getSep;

% GET COMPUTER FOR SPECIAL CASES
Hostname = HF_getHostname;

% Default directories:
DDirs.Root= LF_sysPath;
switch lower(Hostname)
  case 'plethora'; 
    DDirs.Ferrets = 'D:\Data\';
    DDirs.Dropbox = 'C:\Users\dzb\Documents\My Dropbox\';
  case 'genone';
    DDirs.Dropbox = LF_sysPath('Dropbox');
  otherwise
    DDirs.Ferrets = LF_sysPath('project','Records','Ferrets');
    DDirs.Dropbox = LF_sysPath('Dropbox');
end
DDirs.Speaker= LF_sysPath('project','Speakers');
DDirs.G = LF_sysPath('project','Records','G');
DDirs.Kv = LF_sysPath('project','Records','Kv');
DDirs.W = LF_sysPath('project','Records','W');
DDirs.S = LF_sysPath('project','Records','S');
DDirs.SD = LF_sysPath('project','Records','SD');
DDirs.R = LF_sysPath('project','Records','R');
DDirs.M = LF_sysPath('project','Records','M');
DDirs.CBL = LF_sysPath('project','Records','CBL');
DDirs.CBA = LF_sysPath('project','Records','CBA');
DDirs.STRF = LF_sysPath('project','STRFs');
DDirs.WAR = LF_sysPath('project','WARs');
DDirs.WAV = LF_sysPath('project','STRFs','WAV');
DDirs.Records = LF_sysPath('project','Records');
DDirs.bcluster = LF_sysPath('project','Records','bcluster');
DDirs.Export = LF_sysPath('project','Records','Export');
DDirs.Analysa = LF_sysPath('project','Records','Analysa');
DDirs.Notes = LF_sysPath('project','Records','Notes');
DDirs.TIT =  LF_sysPath('project','STRFs','TITemplates');
DDirs.TIIT = LF_sysPath('project','STRFs','TIITemplates');
DDirs.Results = LF_sysPath('project','Records','Results');
DDirs.TmpResults = LF_sysPath('tmp','Results');
DDirs.invitro = LF_sysPath('project','ExtraEPSC');
DDirs.TTNfigures = LF_sysPath('project','TTNdata','Paper','figures');
if ispc DDirs.mfiles = LF_sysPath('Bernhard','mfiles'); 
else DDirs.mfiles = LF_sysPath('mfiles'); end
DDirs.SpikeTimes = LF_sysPath('project','Records','SpikeTimes');
DDirs.VarDelay = LF_sysPath('project','ART_VarDelay');
DDirs.MNTBmodel = LF_sysPath('project','ART_MNTBmodel');
DDirs.MNTBwavedev = LF_sysPath('project','ART_MNTBwavedev');
DDirs.iPPdetect = LF_sysPath('project','ART_iPPdetect');
DDirs.AVCNwaveforms = LF_sysPath('project','ART_AVCNwaveforms');
DDirs.AVCNcluster = LF_sysPath('project','ART_AVCNcluster');
DDirs.AVCNwavefit = LF_sysPath('project','ART_AVCNwavefit');
DDirs.NeuroFEM = LF_sysPath('project','ART_NeuroFEM');
DDirs.ArrayDrive = LF_sysPath('project','ART_ArrayDrive');
DDirs.Colormaps = LF_sysPath('project','STRFs','Colormaps');
DDirs.NoiseCorr = LF_sysPath('project','ART_NoiseCorr');
DDirs.Shepard = LF_sysPath('Dropbox','Projects','Shepard');
DDirs.ShepardPhys = [DDirs.Shepard,'ART_ShepardPhys',Sep];
DDirs.ShepardBehave = [DDirs.Shepard,'ART_ShepardBehave',Sep];
DDirs.StimSpace = LF_sysPath('Dropbox','Projects','ART_StimSpace');
DDirs.ART_MANTA = LF_sysPath('Dropbox','Projects','ART_MANTA');
DDirs.ITBarrel=LF_sysPath('Dropbox','Projects','ART_ITBarrel');
DDirs.Interneuron = LF_sysPath('Dropbox','Projects','IrregularInterneuron');
DDirs.ISH2012 = LF_sysPath('Dropbox','Wissen','Publications','Conferences','ISH','2012','Submission');
DDirs.Texture = LF_sysPath('Dropbox','Psychophysics','ART_TextureDetection');
DDirNames = fieldnames(DDirs);
DDirString = DDirNames{1};
for i=2:length(DDirNames) DDirString = [DDirString,' | ',DDirNames{i}]; end 

% collect and set the User input directories to Dirs (global)
UDirs = cell2struct(varargin(2:2:end),varargin(1:2:end),2);
UDirNames = fieldnames(UDirs);
for i=1:length(UDirNames)
  UDirs.(UDirNames{i}) = sysPathConcat(UDirs.(UDirNames{i}));
  if ~strmatch(UDirNames{i},DDirNames)
    warning(['PathName not in the list of default-path names: ', DDirString]); 
  end
  Dirs.(UDirNames{i}) = UDirs.(UDirNames{i});
end

% set the remaining paths to the default paths, if they have not been set via User or GV
if ~isempty(Dirs)
  for i=1:length(DDirNames)
    if ~sum(strcmp(DDirNames{i},fieldnames(Dirs)))
      Dirs.(DDirNames{i}) = DDirs.(DDirNames{i});
    end
  end
else
  Dirs = DDirs;
end

function Path = LF_sysPath(varargin)
global Dirs; 
if ~isempty(Dirs) & isfield(Dirs,'Root') Root = Dirs.Root; 
else 
  if isunix cd('~/'); Root = pwd;
  elseif strmatch(architecture,'PCWIN') 
    Name = getenv('computername');
    switch Name
      case 'PLETHORA'; Root = 'D:\';
      otherwise Root = 'C:\';
    end
  end
  Dirs.Root = Root;
end
Path = sysPathConcat(Root,varargin);