function testanimate
cont = true;
plot(0, 0, 'o');
set(gcf, 'DeleteFcn', @stopfcn);
theta = 0;
thetainc = 0.1;
radius = 80;
while cont
    plot(0, 0, 'go');
      % compute x and y of orbiting object
      [x, y] = pol2cart(theta, radius);
      hold on;
      plot(x, y, 'ro');
      hold off;
      axis equal;
      axis([-100 100 -100 100]);
      drawnow;
      theta = theta + thetainc;
  end
      function stopfcn(~,~,~)
          cont = false;
      end
  end