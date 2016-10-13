function varargout = undock(handles)
%  VALIDHANDLES = UNDOCK(HANDLES)   undocks the figure(s) with the specified handle(s)
%
% If an output is requested it returns the figureHandlesOrString variable
% (containing only the elements of the original figureHandlesOrString
% variable which are figure handles). The 'windowStyle' property is set to
% 'normal' 
%
%  Examples:
%     
%       undock(aFigureHandle)
%       undock(gcf)
%       undock
%       undock all
%       undock('all')
%       undock(figure(1))
%       undock([figure(1) figure(2) figure(3)])
%   a = undock([figure(1) figure(2) figure(3)])
%
%
%  See also DOCK SETFIGDOCKGROUP GTSCROLL
%
% Written by Riccardo Meldolesi 21-August-2008 


validHandles = dock(handles,[],'undock');

if nargout > 0
   varargout{1} = validHandles;
end




