function HF_matchAspectRatio(varargin)

P = parsePairs(varargin);
if ~isfield(P,'Area') P.Area = 'print'; end
if ~isfield(P,'Factor') P.Factor = 1; end

PP = get(gcf,'PaperPosition');
AR = PP(4)/PP(3);
FP = get(gcf,'Position');

PPI = get(0,'screenpixelsperinch');

% The cmConvFactor is the pixels per cm
% switch lower(HF_getHostname)
%   case {'thot'}; cmConvFactor = 44.9; % assuming MacBook13"    
%   case {'ganesha'}; cmConvFactor = 57; % assuming rMBP 15"
%   case {'genone','deeppurple'}; cmConvFactor = 50.9; % assuming Thinkpad x300
%   otherwise cmConvFactor = 40; warning('Probably need to add correct scaling factor for computer');
% end
cmConvFactor = 44.9;

inConvFactor = 2.54*cmConvFactor ;

switch P.Area
  case 'keep';
    FArea = FP(3)*FP(4);
    Xnew = round(sqrt(FArea/AR));
    Ynew = round(Xnew*AR);
    FP(3:4) = [Xnew,Ynew];
    
  case 'print';
    cUnits = get(gcf,'PaperUnits');
    switch cUnits
      case 'inches'; ConvFactor = inConvFactor;
      case 'centimeters'; ConvFactor = cmConvFactor; 
      otherwise error('Not Implemented for selected Units');
    end
    FP(3:4) = round(ConvFactor*PP(3:4)/P.Factor);
  otherwise
end
     
set(gcf,'Position',FP)
