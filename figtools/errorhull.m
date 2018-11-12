function h = errorhull(X,M,S,varargin)
% Errorhull plots an area around a curve, 
% which is determined by the local standarddeviations at each point
% It thus provides a nices way of diaplaying variability than the command errorbar.
% Arguments are similar to errorbar:
% X : X values
% M : mean of Y values (as Column Vector)
% S : errors in Y direction (as Column Vector, can have two columns for different upper/lower bounds) 
% PlotOpt: linespec style line options 
% HullCol: Color of the error hull (area)
%
% example: errorhull([1:10],rand(1,10)+1,rand(1,10)/5+.5,'r.-',[.8,.8,.8])

P =  parsePairs(varargin);
checkField(P,'Color',[0,0,1]);
checkField(P,'FaceColor',HF_whiten(P.Color,0.5))
checkField(P,'Alpha',0.5);
checkField(P,'LineWidth',1);
checkField(P,'LineStyle','-');

% ONLY OPENGL CAN RENDER TRANSPARENCY
if P.Alpha<1 set(gcf,'Renderer','OpenGL'); end

% CHECK IF UPPER AND LOWER BOUNDS ARE THE SAME OR DIFFERENT
if isempty(find(size(M)==1)) error('Means have to be a vector!'); end
if size(X,1)<size(X,2) X = X'; end; % MAKE COLUMN VECTOR
if size(M,1)<size(M,2) M = M'; end % MAKE COLUMN VECTOR
Dim = find(size(S)==length(M)); % 
if length(Dim)==1 % BOTH DIMENSIONS MATCH
  if Dim ~= 1 S = S'; end
end

if size(S,2)==1
  S(:,2)=S(:,1);
end

NaNInd = union(find(isnan(M)==1),find(isnan(S)==1));
NonNaNInd = setdiff([1:length(M)],NaNInd);
M = M(NonNaNInd); S = S(NonNaNInd,:); X = X(NonNaNInd);
h(2) = curvesurf(X,M-S(:,1),X,M+S(:,2),P.FaceColor); hold on;
set(h(2),'FaceAlpha',P.Alpha,'EdgeColor','none')
if ~strcmp(P.LineStyle,'none')
  h(1) = plot(X,M,'Color',P.Color,'LineWidth',P.LineWidth,'LineStyle',P.LineStyle);
else h(1) = NaN;
end