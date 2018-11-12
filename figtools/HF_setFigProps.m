function HF_setFigProps(varargin)
  
P = parsePairs(varargin);
checkField(P,'FIG',gcf);
checkField(P,'AxisOpt','');
checkField(P,'AxisLabelOpt','');
checkField(P,'TitleOpt','');
checkField(P,'LabelOpt','');

if evalin('caller','exist(''AxisOpt'',''var'')') P.AxisOpt = evalin('caller','AxisOpt'); end; 
if evalin('caller','exist(''LineOpt'',''var'')') P.LineOpt = evalin('caller','LineOpt'); end; 
if evalin('caller','exist(''AxisLabelOpt'',''var'')'); P.AxisLabelOpt = evalin('caller','AxisLabelOpt'); end; 
if evalin('caller','exist(''TitleOpt'',''var'')'); P.TitleOpt = evalin('caller','TitleOpt'); end; 
if evalin('caller','exist(''LabelOpt'',''var'')'); P.LabelOpt = evalin('caller','LabelOpt'); end; 

% ASSIGN AXES PROPERTIES
if ~isempty(P.AxisOpt)
  cAxes = get(P.FIG,'Children');
  Types = get(cAxes,'Type');
  Ind = strcmp(Types,'axes');
  cAxes = cAxes(Ind);
  set(cAxes,P.AxisOpt{:});
end

% ASSIGN AXIS LABEL PROPERTIES
if exist('cAxes','var') && ~isempty(P.AxisLabelOpt)
  for iA=1:length(cAxes)
    cAxisLabelHandles = [get(cAxes(iA),'XLabel'),get(cAxes(iA),'YLabel'),get(cAxes(iA),'ZLabel')];
    set(cAxisLabelHandles,P.AxisLabelOpt{:});
  end
end
  
% ASSIGN AXIS LABEL PROPERTIES
if exist('cAxes','var') && ~isempty(P.TitleOpt)
  for iA=1:length(cAxes)
    cTitleHandles = get(cAxes(iA),'Title');
    set(cTitleHandles,P.TitleOpt{:} );
  end
end
  
% ASSIGN AXIS LABEL PROPERTIES
if exist('cAxes','var') && ~isempty(P.LabelOpt)
  for iA=1:length(cAxes)
    cChildren = get(cAxes(iA),'Children');
    cTypes = get(cChildren,'Tag');
    Ind = strcmp(cTypes,'Label');
    cChildren = cChildren(Ind);
    set(cChildren,P.LabelOpt{:} );
  end
end
  



