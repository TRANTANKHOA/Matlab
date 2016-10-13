function p = time_analysis(myscript)
%time_analysis Profiling the execution time of a given script
%
% _________________________________________________________________________
% INPUT ARGUMENTS
% myscript = handle to the targeted script
% _________________________________________________________________________
% OUTPUT ARGUMENTS
% p        = profile info
%
% -------------------- Copyright (C) 2016 Khoa T. Tran --------------------
%   Created:        05-Aug-2016
%   Last edited:    05-Aug-2016
%   Matlab version: 8.6.0.267246 (R2015b)
%   Email:          khoa.tran@unsw.edu.au
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or (at
%   your option) any later version.
%       This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%       A copy of the GNU General Public License can be found at:
%   http://www.gnu.org/licenses/.
% _________________________________________________________________________

profile on -history
myscript() 
% Save the results as HTML. By default, profsave saves the files to the
% profile_results subfolder in your current working folder.
profsave
% Profile Function and Display Function Call History
% Save the profile results.
p = profile('info');
% Display function entry and exit information by iterating over the
% function call history.
numEvents = size(p.FunctionHistory,2);
for n = 1:numEvents
    name = p.FunctionTable(p.FunctionHistory(2,n)).FunctionName;
    
    if p.FunctionHistory(1,n) == 0
        disp(['Entered ' name]);
    else
        disp(['Exited ' name]);
    end
end