function [  ] = err_handler( err, file )
%ERR_HANDLER Function to save error info out for MCGC
%   Simple Function to save error messages out to a file. 
%   For use with Command Line invocations of MATLAB to make messages
%   persist
%
% USAGE
%
%   err_handler(err,file)
%
% OUTPUT
%
%   Saves the content of an Mexception object err
%   to the file given by file.

%   Jonathan Williams 23/03/2017
%		jonvw28@gmail.com	

fid = fopen(file,'a+');
fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'));
fclose(fid);
end