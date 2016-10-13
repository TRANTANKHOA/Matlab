%# ellipse centered at (0,0) with axes length
%# major=20, ,minor=10, rotated 50 degrees
%# (drawn using the default N=36 points)
clf
axis normal
hold on
p = calculateEllipse(0, 0, 1.6, 16, 60);
fill(p(:,1), p(:,2),[0.9 0.9 0.9])
p = calculateEllipse(0, 0, 0.8, 8, 60);
fill(p(:,1), p(:,2),[0.8 0.8 0.8])
p = calculateEllipse(0, 0, 0.2, 2, 60);
fill(p(:,1), p(:,2),[0.7 0.7 0.7])
p = calculateEllipse(0, 0, 0.05, 0.5, 60);
fill(p(:,1), p(:,2),[0.6 0.6 0.6])

%% Draw slices
grid on
line([5 5],[0 6])
line([4.7 5.3],[2 2])
line([4.7 5.3],[4 4])
line([2 10],[3.5 3.5])
line([4 4],[3.2 3.8])
line([6 6],[3.2 3.8])
line([8 8],[3.2 3.8])
xlabel('$$x_1$$')
ylabel('$$x_2$$')
text(11,2.5,'Slice 1')
text(1,7.5,'Slice 2')

%% Draw eige-axes
grid on
line([-15 15],tan(pi/6)*15*[-1 1],'LineStyle','--')
line([-1.5 1.5],tan(-pi/3)*1.5*[-1 1],'LineStyle','--')
line([-0.5 1.5]+10,tan(-pi/3)*1*[-1 1]+10*tan(pi/6))
line([-14 16],tan(pi/6)*15*[-1.01 0.99])
text(12,4,'Slide 1')
text(17,9,'Slide 2')
xlabel('$$x_1$$')
ylabel('$$x_2$$')
text(-10,8,'Principal axes')
grid off







%% Draw boxes
grid on
line([-8 -8],[-4 -6])
line([-8 -6],[-4 -4])
line([-8 -6],[-6 -6])
line([-6 -6],[-4 -6])
line([-8 -8],[-8 -6],'LineStyle',':')
line([-8 -8],[-4 -2],'LineStyle',':')
line([-6 -6],[-4 -2],'LineStyle',':')
line([-8 -10],[-4 -4],'LineStyle',':')
line([-6 -4],[-4 -4],'LineStyle',':')
line([-8 -10],[-6 -6],'LineStyle',':')
line([-10 -10],[-2 -8])
line([-4 -4],[-2 -8])
line([-4 -10],[-2 -2])
line([-4 -10],[-8 -8])
line([-4 -2],[-2 -2],'LineStyle',':')
line([-4 -4],[0 -2],'LineStyle',':')
line([-2 0],[0 0],'LineStyle',':')
line([-2 -2],[0 2],'LineStyle',':')
line([-10 -10],[-2 0])
line([-2 -10],[0 0])
line([-2 -2],[0 -8])
line([-2 -4],[-8 -8])
line([-2 0],[-8 -8])
line([0 0],[-8 2])
line([-10 -10],[0 2])
line([-10 0],[2 2])
text(0.5,-7,'Multivariate slice')
grid off
