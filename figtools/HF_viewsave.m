function HF_viewsave(varargin)
% Helper function for saving and viewing figure output
%  
% Options:
% -name : file name (without extension)
% -path : file path (if left out, set to empty)
% -view [{0},1] :
% -format [{epsc2},epsc,pdf,jpg,tif,...] : graphics file format 
% -tag [{=format}] : file extension 
% -res [{400}] : dpi resolution
% -fullpage [{0},1] : places the figure on a full A4 page
% -label [String] : places a figure label underneath the figure

% Dirs = setgetDirs; 
if ispc Sep = '\'; else Sep = '/'; end

% PARSE ARGUMENTS
if isstr(varargin{3}) % new format 
  P = parsePairs(varargin);
else % old format
  P.name = varargin{2}; P.view = varargin{3};
  P.path = [Dirs.Root,'project',Sep,varargin{1},Sep,'figures'];
  if length(varargin)>1 P.format = varargin{4}; end
  if length(varargin)>2 P.tag = varargin{5}; end
  if length(varargin)>3 P.res = varargin{6}; end
end

if P.save

% SET DEFAULT ARGUMENTS
checkField(P,'name');
checkField(P,'view',0);
checkField(P,'format','pdf');
checkField(P,'tag',P.format);
checkField(P,'res',400);
checkField(P,'path',[]);
checkField(P,'fullpage',0);
checkField(P,'label',[]);
checkField(P,'colorspace','');
if ~isempty(P.path) & ~P.path(end)==Sep    P.path(end+1) = Sep; end
outfile = [P.path,P.name];

% PLACE FIGURE ON A FULL PAGE
global Fullpage; if Fullpage P.fullpage = 1; end
if P.fullpage
  if length(P.format)>2 & P.format(1:3)=='eps' 
    P.format = P.format(2:end); % avoid EPS, since it would not allow full pages
  end
  set(gcf,'PaperType','A4'); PS = get(gcf,'PaperSize');
  PP = get(gcf,'PaperPosition');
  PP = [(PS-PP(3:4))/2,PP(3:4)]; set(gcf,'PaperPosition',PP);
end

% ADD FIGURE LABEL
% Requires some Postscript modifications to work
% first search for (label) e.g. (Figure 1) s
% get Position (stated above it), e.g. 1506 2369 mt 
% move down (increase second value)
% and change viewing box at beginning of PS, e.g.   291   127  2643  2380 rc
% In general the position format is xstart,ystart,xstop,ystop
% Probably implement this in a Perl script
if ~isempty(P.label)
  PP = get(gcf,'PaperPosition');
  if PP(2)<0.5 % MOVE FIGURE UPWARDS
    PP(2) = 0.5; set(gcf,'PaperPosition',PP);
  end
  % AP = [PP(1)+PP(3)/2,PP(2)-0.5,PP(3),0.4];
  AP = [0.5,-.5,.1,.1];
  h=axes('OuterPosition',[0,0,1,1],'visible','off');
  text(AP(1),AP(2),P.label,'Units','normalized',...
    'FontSize',14,'FontWeight','bold',...
    'HorizontalAlignment','Center');
end

% SET SPECIAL OPTIONS DEPENDING ON THE RENDERER
switch lower(get(gcf,'Renderer'))
  case 'painters'; 
  otherwise if ~isfield(P,'res') | isempty(P.res) P.res = 400; end
end

% SET TAG FOR FILETYPE
P.format = lower(P.format);
switch P.format
  case {'ps','psc','psc2'};    P.tag = 'ps';
  case {'eps','epsc2','epsc'}; P.tag = 'eps';
  case {'meta'}; P.tag = 'emf';
  otherwise P.tag = P.format;
end
outfiletag = [outfile,'.',P.tag];

% DELETE PREVIOUS VERSION OF OUTPUT FILE
if exist(outfiletag,'file') delete(outfiletag); end

% PRINT THE FIGURE
%'-loose' : could be useful for EPS printing 
print(gcf,['-d',P.format],[P.colorspace],['-r',num2str(P.res)],[outfile,'.',P.tag]); 

% CHECK FOR VIEWER PROGRAMS
switch computer
  case {'MAC','MACI','MACI64'};
  case 'PCWIN';
    isDistiller = exist('C:\Programme\acroread.lnk','file');
    isReader = exist('C:\Programme\acroread.lnk','file');
end

% CONVERT TO PDF IF FILE IS EPS
if strcmp(P.tag,'eps')
  if isunix
    if exist([outfile,'.pdf'],'file') delete([outfile,'.pdf']); end
    system(['pstopdf ',outfiletag,' ',outfile,'.pdf']);
  else
    if isDistiller  system(['C:\Programme\acrodist.lnk ',outfiletag]);
    else LF_DistError;
    end
  end
end

% VIEW FILE
if P.view
  switch computer
    case {'MAC','MACI','MACI64'};
        outfile=strrep(outfile,' ','\ ');
      switch P.format
        case {'ps','psc','psc2','eps','epsc','epsc2','pdf'}
          system(['osascript -e "tell application \"Preview\" to close document \"',P.name,'.pdf\" " 2>/dev/null ']);
          system(['osascript -e "tell application \"Preview\" to close document \"',outfile,'.pdf\" " 2>/dev/null ']);
          system(['open -a Preview ',outfile,'.pdf']);
        otherwise
          system(['osascript -e "tell application \"Preview\" to close document \"',P.name,'.',P.tag,'\" " 2>/dev/null']);
          system(['osascript -e "tell application \"Preview\" to close document \"',outfile,'.',P.tag,'\" " 2>/dev/null ']);
          system(['open -a Preview ',outfile,'.',P.tag]);
      end
      
    case 'PCWIN';
      switch P.format
        case {'ps','psc','psc2','eps','epsc','epsc2','pdf'}
          if isReader  system(['C:\Programme\acroread.lnk ',outfile,'.pdf']);
          else LF_ReaderError
          end
        otherwise
          fprintf('No Viewer Programm set!');
      end
      
    case {'GLNX','GLNXA64'}
      switch P.format
        case 'pdf'
          [s,e] = system(['evince -w ',outfile,'.pdf']);
        otherwise
          [s,e] = system(['eog -g ',outfile,'.',P.format,' &']);
      end
      
    otherwise fprintf(['Computertype not implemented: ',computer,'\n']);
  end
end

end
function LF_DistError
  error('EPS file cannot be converted to PDF, since Distiller was not found.');
  
function LF_ReaderError
  error('PDF cannot be displayed, since the Acrobat Reader was not found.');