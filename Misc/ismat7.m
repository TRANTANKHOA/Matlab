function tf = ismat7(fname)
% Check if the version of the Matlab data file as 7.3 MAT-file
	x = evalc(['type(''', fname, ''')']);
	tf = strcmp(x(2:20), 'MATLAB 7.3 MAT-file');
end