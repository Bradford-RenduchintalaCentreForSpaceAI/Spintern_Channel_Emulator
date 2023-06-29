close all 
clearvars
g = hgtransform;
q = animatedline;
x = [0 1 1 0];
y = [0 0 1 1];
patch('XData',x,'YData',y,'FaceColor','red','Parent',g)
axis equal
xlim([-10 10])
ylim([-10 10])


pt1 = [-10 0 0];
pt2 = [9  0 0];
for t=linspace(0,1,100)
  g.Matrix = makehgtform('translate',pt1 +t*(pt2-pt1));
  drawnow
end

% Speed of reciever stays still 
% At freq of 1MHz
figure

x1 = linspace(0,100,10000);


plot(x1,sin(x1*(3E8/(3E8))))
