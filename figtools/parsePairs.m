function S = parsePairs(Cell)

%Structures varargin input and checks for invalid input. The input has to
%come in 'argument','value' pairs, where 'argument' is a variable of type 
%char. The parameters are extracted from the varargin cell and ordered in
%the structure S.
S = [];
if mod(length(Cell),2)==1 error('Number of input arguments is not even!'); end
for i=1:length(Cell)/2
 if ~ischar(Cell{2*i-1}) error('One of the cells is not a Name.'); end
 eval(['S.',Cell{2*i-1},' = Cell{2*i};']);
end