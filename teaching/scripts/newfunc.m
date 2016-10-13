function newfunc(name)
%NEWFUNC Customize and create a formatted file for a new MATLAB function
%
%  USAGE:   newfunc(name)
%           *newfunc with no args prints help and shows default VARARGINs
% 
% _________________________________________________________________________
%  NECESSARY ARGUMENT
%     name          = function name e.g. ['myfunc.m','myfunc',{'myfunc'}] 
% _________________________________________________________________________
%  EXAMPLES
%     newfunc
%     newfunc('mynewfunc'); 
% 
% ---------------------- Copyright (C) 2014 Bob Spunt ----------------------
%	Created:  2014-09-27
%	Email:    spunt@caltech.edu
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or (at
%   your option) any later version.
%       This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%       You should have received a copy of the GNU General Public License
%   along with this program.  If not, see: http://www.gnu.org/licenses/.
% _________________________________________________________________________

%% Printing help or check for duplication
if nargin==0, help newfunc; return; end
name = [name '.m'];
if exist(name, 'file')
    resp = inputoption(['File by that name already exists at: \n'...
        which(name) '\nDo you want to proceed?'], {'y' 'n'}); 
    if strcmp(resp, 'n'), return; end
end
%% CREATE OPEN-MFILE
[~, path] = uiputfile(name,'Save as');
path      = fullfile(path,name);
fid       = fopen(path,'w+');
%% Template
today = date; thisyear = year(today);
fprintf(fid, ['function [] = ' name(1:end-2) '()\n']);
fprintf(fid, '%s\n' ...
    ,['%' name(1:end-2) ' ']...
    ,'%'...
    ,'% _________________________________________________________________________'...
    ,'% INPUT ARGUMENTS'...
    ,'%'... 
    ,'% _________________________________________________________________________'...
    ,'% OUTPUT ARGUMENTS'...
    ,'%'... 
   ,['% -------------------- Copyright (C) ' num2str(thisyear) ' Khoa T. Tran --------------------']...
   ,['%   Created:        ' today]...
   ,['%   Last edited:    ' today]...
   ,['%   Matlab version: ' version]...
    ,'%	  Email:          khoa.tran@unsw.edu.au'...
    ,'%'...
    ,'%   This program is free software: you can redistribute it and/or modify'...
    ,'%   it under the terms of the GNU General Public License as published by'...
    ,'%   the Free Software Foundation, either version 3 of the License, or (at'...
    ,'%   your option) any later version.'...
    ,'%       This program is distributed in the hope that it will be useful, but'...
    ,'%   WITHOUT ANY WARRANTY; without even the implied warranty of'...
    ,'%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU'...
    ,'%   General Public License for more details.'...
    ,'%       A copy of the GNU General Public License can be found at:'...
    ,'%   http://www.gnu.org/licenses/.'...
    ,'% _________________________________________________________________________'...
    );

%% Done
fclose(fid);
fprintf('\nFunction saved to: %s\n', name); 
edit(path);
end
 
function resp       = inputoption(prompt, respset)
% INPUTOPTION Like INPUT, but with default and response option set
%
%  USAGE: resp = inputoption(prompt, respset)
% __________________________________________________________________________
%  INPUTS
%	prompt:     message to user 
%	respset:    valid responses (case is ignored when validating)
%
% __________________________________________________________________________
%  EXAMPLES
%	>> resp = inputoption('Would you like to quit?', {'Y' 'N'}, 'Y')
%	>> resp = inputoption('1, 2, or 3 apples?', [1 2 3])
%

% ---------------------- Copyright (C) 2014 Bob Spunt ----------------------
%	Created:  2014-09-29
%	Email:    spunt@caltech.edu
% __________________________________________________________________________
if nargin < 2, disp('USAGE: resp = inputoption(prompt, respset)'); return; end
if iscell(prompt), prompt = char(prompt); end
numtag = 0; 
if any(isnumeric(respset))
    numtag = 1; 
    if size(respset, 1)==1, respset = respset'; end
    respset = cellstr(num2str(respset)); 
end
if ischar(respset), respset = cellstr(respset); end
nopt = length(respset);
respset = strtrim(respset); 
if nopt > 2
    optstr = sprintf(['(' repmat('%s, ', 1, nopt-1) 'or %s)'], respset{:});
else
    optstr = sprintf('(%s or %s)', respset{:}); 
end
resp = input(sprintf('%s %s ', prompt, optstr), 's'); 
while ~any(strcmpi(respset, resp))
    resp = input(sprintf('Invalid response. %s ', optstr), 's'); 
end
if numtag, resp = str2double(resp); end
end
