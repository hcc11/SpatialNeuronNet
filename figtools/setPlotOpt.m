function setPlotOpt(Project,varargin)

P = parsePairs(varargin);
% Make General Definitions (can be overridden in the following)
FontName = 'Arial'; LocalUnits = 'centimeters';
gvec = [1,1,1]; if ispc Sep = '\'; else Sep = '/'; end
MaxHeight = 23.5; % A4 MaxHeight
inpathS = '[P.path,''figures'',Sep]'; 
outpathS = '[P.path,''Submission_'',Project,Sep]';
Renderer = 'Painters'; PaperPositionMode = 'Manual';
if ~isfield(P,'width') P.width = 4; end; if ~isfield(P,'height') P.height = 3; end

switch lower(Project)
  case 'custom' % Define things by Hand
    LabelFont = {12}; TitleFont = {12}; AxisLabelFont = {12}; 
    AxisFont = {8}; LegendFont = {8}; TextFont = {10};   
    Cols = [1,1.5,2;8.25,12.5,17.15];
    AxisOpt = {'LineWidth',1};
%     inpathS = '[P.inpath,Sep]';
%     outpathS = '[P.outpath,Sep]';
    
    % BY PRESENTATION TYPE
  case 'manuscript'
    LabelFont = {12}; TitleFont = {9}; AxisLabelFont = {8};
    AxisFont = {7}; LegendFont = {8}; TextFont = {6};
    Cols = [1,1.5,2;8.25,12.5,17.5]; PointOpt = {'MarkerSize',8};
    AxisOpt = {'LineWidth',1};LineOpt = {'LineWidth',1};
%     outpathS = '[P.path,''Manuscript'',Sep]';

  case 'a4paper'
    LabelFont = {10}; TitleFont = {9}; AxisLabelFont = {8}; AxisFont = {6}; LegendFont = {8};
    outpathS = '[P.path,''figures'',Sep]'; 
    
  case 'poster'
    LabelFont = {18}; TitleFont = {22}; AxisLabelFont = {18}; AxisFont = {14}; LegendFont = {16};
    AxisOpt = {'LineWidth',1.5}; TextFont = {14}; 
    outpathS = '[P.path,''figures'',Sep]'; clear MaxHeight;

  case 'presentation'
    LabelFont = {12}; TitleFont = {20}; AxisLabelFont = {14}; AxisFont = {12}; LegendFont = {14};    
    outpathS = '[P.path,''figures'',Sep]'; TextFont = {12}; AxisOpt = {'LineWidth',1.5};
    
  % BY PROJECT
  case 'mntbmodel'
    LabelFont = {10}; TitleFont = {8}; AxisLabelFont = {7}; AxisFont = {6}; LegendFont = {9};    
    
  % BY JOURNAL
  case 'jneurosci'
    LabelFont = {12}; TitleFont = {10}; AxisLabelFont = {9}; AxisFont = {8}; LegendFont = {9}; 
    Cols = [1,1.5,2;8.5,12.5,17.15];

  case 'neuron'
    LabelFont = {11,'bold'}; TitleFont = {10}; AxisLabelFont = {8}; AxisFont = {7}; LegendFont = {8};
    TextFont = {7}; AxisOpt = {'LineWidth',1};
    Cols = [1,1.5,2;8.5,12.5,17.15]; 
  
  case 'pnas'
    LabelFont = {12}; TitleFont = {10}; AxisLabelFont = {9}; AxisFont = {8}; LegendFont = {9};
    TextFont = {8};
    Cols = [1,1.5,2;8.5,12.5,17.15]; 
  
  case 'jaro'
    LabelFont = {12}; TitleFont = {10}; AxisLabelFont = {9}; AxisFont = {8}; LegendFont = {9};    
    Cols = [1,1.5,2;8.5,12.5,17.15];

  case 'jnphys'
    LabelFont = {12}; TitleFont = {10}; AxisLabelFont = {9}; AxisFont = {8}; LegendFont = {9};    
    Cols = [1,1.5,2;8.5,12.7,18];

  case 'plos'
    LabelFont = {12}; TitleFont = {10}; AxisLabelFont = {9}; AxisFont = {8}; LegendFont = {9};    
    Cols = [1,1.5,2;8.25,12.5,17.15];
    
  case 'jasa'
        LabelFont = {12}; TitleFont = {10}; AxisLabelFont = {9}; AxisFont = {8}; LegendFont = {9};
        Cols = [1,1.5,2;8.5,13,17.15];
       
  case 'finc'
    LabelFont = {12}; TitleFont = {10}; AxisLabelFont = {9}; AxisFont = {8}; LegendFont = {9};    
    Cols = [1,1.5,2;8.25,12.5,17.15];

  case 'nc'
    LabelFont = {10}; TitleFont = {9}; AxisLabelFont = {8}; AxisFont = {7}; LegendFont = {8};
    Cols = [1,1.5,2;10.8,NaN,NaN];

  % BY CONFERENCE
  case 'aro'
    LabelFont = {22}; TitleFont = {20}; AxisLabelFont = {18}; AxisFont = {14}; LegendFont = {16};    
    AxisOpt = {'LineWidth',1.5}; LineOpt = {'LineWidth',2}; MarkerOpt = {'MarkerSize',20};
    AxesOpt = {};
    outpathS = '[P.path,''ARO_Poster'',Sep]';

  case 'ish';
    LabelFont = {12}; TitleFont = {10}; AxisLabelFont = {7}; AxisFont = {6}; LegendFont = {7}; TextFont = {6};
    AxisOpt = {'LineWidth',0.5}; LineOpt = {'LineWidth',1}; MarkerOpt = {'MarkerSize',20};
    AxesOpt = {};
    outpathS = inpathS;    
    
  otherwise error('Project not defined');
end

% Perform automatic assignments
inpath = eval(inpathS); outpath = eval(outpathS);
if ~exist('FigureSize','var') & isfield(P,'height') & isfield(P,'cols') & exist('Cols','var')
  FigureSize = [Cols(2,find(P.cols==Cols(1,:))),P.height];
else FigureSize = [P.width,P.height]; end
if ~exist('PaperSize','var') PaperSize = FigureSize; end
if ~exist('PaperPosition','var') PaperPosition = [0,0,FigureSize]; end
if isfield(P,'Renderer') Renderer = P.Renderer; end
FigOpt = {'PaperPositionMode',PaperPositionMode,...
  'PaperUnits',LocalUnits,'PaperSize',PaperSize,'PaperPosition',PaperPosition,...
  'Renderer',Renderer,'Color','white'};

FontVars = whos('*Font');
for i=1:length(FontVars)
  cVar = eval(FontVars(i).name); VarName = [FontVars(i).name(1:end-4),'Opt'];
  if exist(VarName,'var') tmp = eval(VarName); 
  elseif isfield(P,VarName) tmp = P.(VarName);
  else tmp = {}; end
  if ~isempty(cVar{1}) tmp = ['FontSize', cVar{1}, tmp];  end
  if length(cVar)>=2 & ~isempty(cVar{2}) tmp = ['FontWeight', cVar{2}, tmp];  end
  if length(cVar)>=3 & ~isempty(cVar{3}) tmp = ['FontName', cVar{3}, tmp];  
  else tmp = ['FontName', FontName, tmp];
  end; eval([VarName,'=tmp;']);
end

% Test options
if exist('MaxHeight','var') LF_CheckHeight(P.height,MaxHeight,LocalUnits); end
% if ~exist(inpath,'file') error(['Path ',inpath,' does not exist!']); end
% if ~exist(outpath,'file') error(['Path ',outpath,' does not exist!']); end

% Assign variables in caller Workspace
Vars = {'gvec','Sep','inpath','outpath',...
  'FigOpt','TitleOpt','LabelOpt','AxisLabelOpt','AxisOpt','TextOpt',...
  'LegendOpt','LineOpt','MarkerOpt','PointOpt'};
for i=1:length(Vars)
  if exist(Vars{i},'var')  eval(['assignin(''caller'',Vars{i},',Vars{i},');']);  end
end
if isfield(P,'FIG') 
  figure(P.FIG); clf; set(gcf,FigOpt{:});
  HF_matchAspectRatio;
end
if isfield(P,'Bare') & P.Bare == 1; 
  set(gcf,'MenuBar','none','Toolbar','none');
end
Stack = dbstack; assignin('caller','name',Stack(2).name);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = LF_CheckHeight(Height,Max,Units)
if Height>Max warning(['Height must be smaller than ',n2s(Max),' ',Units]); end

function W = LF_Cols(Cols,SelCol)
ind = find(Cols(1,:)==SelCol);
if ~isempty(ind)  W = Cols(2,ind);
else error('Desired Columns not available!');
end