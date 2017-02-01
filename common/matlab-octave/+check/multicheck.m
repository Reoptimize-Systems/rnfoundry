function multicheck (checkfcn, errormsg, errid, varargin)
% performs check on multiple inputs and raises error on failure
%
% Syntax
%
% multicheck (checkfcn, errormsg, errid, arg1, arg2, ..., argn )
%
% Input
%
%  checkfcn - function to apply to arguents to check their validity
%
%  errormsg - error message to display in the event of any failure to pass
%    the check.
%
%  errid - error id to use for failure error message
%
%  arg1, arg2, ..., argn - any number of arguments for which to perform the
%    test.
%

    result = cellfun (checkfcn, varargin, 'UniformOutput', true);
    
    if any (result == false)
        failedinds = find (result ==  false);
        
        failmessage = sprintf ('inputs %s failed validation with following error message:\n%s', ...
                    mat2str (failedinds), errormsg );
                 
        if isempty (errid)
            error (failmessage);
        else
            error (errid, failmessage);
        end
    end
    

end