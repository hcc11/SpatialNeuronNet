function h=FigLabel(string,pos,varargin)

AxisPos = get(gca,'Position');
pos = [pos(1)/AxisPos(3),pos(2)/AxisPos(4)];
pos(2) = pos(2)+1;

h=text(0,0,string); 
set(h,'Units','normalized','HorizontalAlignment','Left',...
  'Position',pos,varargin{:},'Tag','Label');
