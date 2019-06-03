    close all;
    clear all;
    clc;
    theta=0; w=20/60; t=0:10:1440;
     c=0;   env=0.00001; 
v = VideoWriter('newfile2.avi', 'Uncompressed AVI');
v.FrameRate = 10;  % Default 30
%v.Quality = 100;  
open(v)
while c<2000
        %Show fast you want it 
        theta=theta+20;
        %y_plot=y*sind(theta);
        x=0:10:1440;
        x=x+theta;
        y=sind(x + w*t);
        i=ceil(random('uniform',10,20,[1]));
        h1=plot(t,y,'k-','LineWidth',2);
                axis([100,1000,-1.2,1.2]);
          set(gcf,'color','w')   
         set(gca,'xtick',[],'YTick',[],'xcolor',[1 1 1],'ycolor',[1 1 1])
         drawnow 
      %pause(0.0001)
       frame = getframe(gcf);
        writeVideo(v,frame);
	     c=c+1;
        %myMovie(c) = thisFrame;
if (mod(c,100)==0)
        for ix = 15:10:floor(random('uniform',t(8),t(12),1))
    hold on
    plot(t(ix),y(ix),'d', 'MarkerSize', 10,'MarkerEdgeColor','black', ...
        'MarkerFaceColor',[1,0.647,0]);    
    drawnow    
    pause(0.07)     
     frame = getframe(gcf);
        writeVideo(v,frame);
    hold off
    %pause(0.2)                                % control animation speed
    %snapnow                                   % (required for published document only)
        end
end
hold off
end
   %hold off

close(v)
	%close(writerObj);