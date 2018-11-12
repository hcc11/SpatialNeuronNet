function [varargout]=raster2D_ani(s0,t0,t1,Ne1)
% generate movie of spike rasters from 2D spatial network 
% s0: spike train data, s(1,:): spike times, s(2,:): neuron ID  
% t0: starting time of the movie 
% t1: end time
% Ne1: # of neurons per dimension 
% varargout is a struct of frames 

    dta=2; % bin size for raster animation
    timea=t0:dta:t1; % Timeframes
%     timea=dta:dta:T; % Timeframes
%     spkduration=1; % Duration of each spike in raster
    fig1=figure;
    set(gcf,'color','w')
    set(gcf,'position',[350 100 550 450])
%     winsize = get(fig1,'Position');
%     winsize(1:2) = [0 0];
    numframes=numel(timea);
%     A=moviein(numframes,fig1,winsize);
    A(1:numframes) = struct('cdata', [],'colormap', []);
    for i=1:numframes

      % Find exc spikes in this time bin  
      if(i==1)
          Is=find(s0(1,:)<=timea(i)& s0(1,:)>timea(i)-dta & s0(2,:)>0);
          if size(s0,1)==2
              x=ceil(s0(2,Is)/Ne1);
              y=mod(s0(2,Is)-1,Ne1)+1;
          else
              x=s0(2,Is);
              y=s0(3,Is);
          end
          if isempty(Is)
              x=[0];y=[0];
          end
          h=plot(x,y,'k.','MarkerSize',5); % plot command
          set(h,'XDataSource','x');
          set(h,'YDataSource','y');
          axis([0 Ne1 0 Ne1])
          xlabel('neuron location (X)','fontsize',18)
          ylabel('neuron location (Y)','fontsize',18)
          set(gca,'fontsize',15)
          set(gca,'xtick',[0 100 200])
          set(gca,'ytick',[0 100 200])
      else
          Is=find(s0(1,:)<=timea(i) & s0(1,:)>timea(i)-dta & s0(2,:)>0);
      end
      
      % Plot Raster
      if size(s0,1)==2
          x=ceil(s0(2,Is)/Ne1);
          y=mod(s0(2,Is)-1,Ne1)+1;
          pause(.1)
      else
          x=s0(2,Is);
          y=s0(3,Is);
      end
      
      refreshdata(h,'caller')
      drawnow
      
      title(sprintf('t=%d msec',round(timea(i))-t0))
      A(i)=getframe(fig1);
    end
    if nargout==1 
        varargout{1}=A;
    end