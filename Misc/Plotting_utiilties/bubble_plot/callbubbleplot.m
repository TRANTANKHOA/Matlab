% script file to exercise Bubbleplot.m
%%
x1 = 10*rand(10,1);   % generate a first dataset
y1 =  5*rand(10,1);
z1 =  7*rand(10,1);
color1 = [0 0 1];    % make plot symbols blue 
sf1 = 25;
legwords1 = 'Dataset #1';
%-------------------------------------------------
x2 = 18*rand(10,1);   % generate a second dataset
y2 = 8*rand(10,1);
z2 =  5*rand(10,1);
color2 = [1 0 0];    % make plot symbols red
sf2 = sf1*max(z2)/max(z1);             % adjust scale factor for second plot for continuity with first plot
legwords2 = 'Dataset #2';
%-------------------------------------------------
x3 = 7*rand(10,1);   % generate a third dataset
y3 = 6*rand(10,1);
z3 =  9*rand(10,1);
color3 = [0 1 1];    % make plot symbols cyan 
sf3 = sf1*max(z3)/max(z1);             % adjust scale factor for third plot for continuity with first & second
legwords3 = 'Dataset #3';
%-------------------------------------------------
close all
line(x1(1),y1(1),'Color',color1,'Visible','off');  % plot only the first point from each dataset so can use
hold on                                            % legend
line(x2(1),y2(1),'Color',color2,'Visible','off');
line(x3(1),y3(1),'Color',color3,'Visible','off');

h1 = bubbleplot(x1,y1,z1,color1,sf1);  % plot first dataset
hold on;                                      % note hold must be turned back 'on' after each call to bubbleplot 
h2 = bubbleplot(x2,y2,z2,color2,sf2);  % plot second dataset    
hold on;
h3 = bubbleplot(x3,y3,z3,color3,sf3);  % plot third dataset

set(gca,'FontSize',8);
xlabel('Variable #1');
ylabel('Variable #2');
title('3-Dataset Bubbleplot Example   :::   Symbol Size is Proportional to Z-variable Magnitude');
legend(legwords1,legwords2,legwords3,'Location','Best');
set(gcf,'Color',[.8 .8 .8],'InvertHardCopy','off');
%%