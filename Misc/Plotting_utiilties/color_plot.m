function [ h ] = color_plot( x, Y)
%This helps you plot columns of matrix Y over vector x and format the ticks
%   Detailed explanation goes here
n_lines=size(Y,2);
ColorSet = varycolor(n_lines);
hold on;
for i=1:n_lines
    plot(x,Y(:,i),'Color', ColorSet(i,:));
end
h=1;

