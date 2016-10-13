% Dock all figure to main window.
global font_size
font_size=14;
% set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultFigureWindowStyle','normal')
set(0,'defaultUIControlFontName', 'Times New Roman')
% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', font_size)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', font_size)

%Use Latex as default typeset
set(0,'DefaultTextInterpreter','Latex') 

%White figure background
set(0,'DefaultFigureColor','w');

%Set default figure size
w = 1;h = 1;
w=round(w*1050);h=round(h*800);
sz=[100 100 w h];
set(0, 'DefaultFigurePosition', sz);

% %Render plot
% set(0,'DefaultFigureRenderer', 'painters');
set(0,'DefaultLineLineWidth',1.1);
% set(0,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder','-|--|:|-.|.')
% %Change current directory
% % cd AUCC2014
% 
% % Turn off the graphic smoothing feature on MATLAB 2014b onward
% set(groot,'DefaultFigureGraphicsSmoothing','off')
