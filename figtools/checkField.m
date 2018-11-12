function checkField(P,Field,DefVal,Type)

persistent AllFields;

if nargin==1 % CHECK FOR FIELDS THAT DO NOT EXIST (ARE NOT TESTED ABOVE)
  WrongFields = setdiff(fieldnames(P),AllFields);
  if ~isempty(WrongFields) error(['Field ',WrongFields{1},' does not exist!']); 
  end
  return;
else
  AllFields{end+1} = Field; AllFields = unique(AllFields);
end

if ~isempty(P) FN = fieldnames(P); else FN = {}; end
ApproxMatch = strcmp(lower(FN),lower(Field));
if  sum(ApproxMatch)% A POTENTIALLY INEXACT MATCH EXiSTS
  if ~sum(strcmp(FN,Field)) % NO EXACT MATCH, CORRECT TO
    UserField = FN{ApproxMatch};
    fprintf(['Argument name ''',UserField,''' corrected to ''',Field,'''\n']);
    P(1).(Field) = P.(UserField);
    P = rmfield(P,UserField);
  end
else % NO MATCH, NOT EVEN INEXACT EXISTS
  if ~exist('DefVal','var')
    error(['checkField : The argument ''',Field,''' needs to be assigned!']);
  else
    P(1).(Field) = DefVal;
  end
end
assignin('caller','P',P);