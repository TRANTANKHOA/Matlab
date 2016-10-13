% FINISHSAV  Save workspace variables  
%   Change the name of this file to FINISH.M 
%   and put it anywhere on your MATLAB path.
%   When you quit MATLAB this file will be executed.
%   This script saves all the variables in the 
%   work space to a MAT-file.  

%   Copyright 1984-2000 The MathWorks, Inc. 
%   $Revision: 1.4.2.1 $  $Date: 2011/01/28 18:50:30 $

disp(getString(message('MATLAB:finishsav:SavingWorkspaceData')));
save 
