function FIGs = placeFigures(varargin)
% Automatically place figures in a grid (underlying HF_axesDivide allows noneven splits)

%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'Figures','all');
checkField(P,'XY','auto');
checkField(P,'Toolbar','none');
checkField(P,'Menubar','none');
checkField(P,'Position','screen');
checkField(P,'XMargin',0.02);
checkField(P,'YMargin','auto');

%% GET FIGURES TO BE ORGANIZED 
FIGs = get(0,'Children');
if ~strcmp(P.Figures,'all') FIGs = P.Figures; %intersect(FIGs,P.Figures); 
else FIGs = sort(FIGs,'ascend');
end

%% ADAPT GRID TO AVAILABLE FIGURES
NFigs = length(FIGs);
if ~iscell(P.XY) 
  if strcmp(P.XY,'auto')
    P.XY = [ceil(sqrt(NFigs)),ceil(sqrt(NFigs))];
  end
  YHeight = P.XY(2);
  NGrid = P.XY(1)*P.XY(2);
else % GRID IS A CELL, ASYMMETRIC DIVISION
  NGrid = length(P.XY{1})*length(P.XY{2});
  YHeight = mean(P.XY{2});
end
if NGrid < NFigs
  fprintf('WARNING : Adapting grid to the number of figures!\n'); 
  P.XY(1) = ceil(NFigs/P.XY(1)); 
end

%% ADAPT MARGINS
Opts = {}; MAdd = 1;
if ~isempty(P.Toolbar) Opts = {Opts{:},'Toolbar',P.Toolbar}; MAdd = MAdd+1; end
if ~isempty(P.Toolbar) Opts = {Opts{:},'Menubar',P.Menubar}; MAdd = MAdd + 1; end

%% CREATE FIGURE POSITIONS
if strcmp(P.YMargin,'auto') P.YMargin = 0.07*MAdd*YHeight; end
if strcmp(P.Position,'screen') P.Position = get(0,'ScreenSize') + [1,1,-2,-80];end
if ~iscell(P.XY)
  Positions = HF_axesDivide(P.XY(1),P.XY(2),P.Position,P.XMargin,P.YMargin);
else
  Positions = HF_axesDivide(P.XY{1},P.XY{2},P.Position,P.XMargin,P.YMargin);
end
Positions = Positions';

%% POSITION FIGURES
for i=1:NFigs
%  figure(FIGs(i)); % raise figures : too slow 
  cOpts = {Opts{:},'Position',round(Positions{i})};
  set(FIGs(i),cOpts{:});
end

if ~nargout clear FIGs; end 