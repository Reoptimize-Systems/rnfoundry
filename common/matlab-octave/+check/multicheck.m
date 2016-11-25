function multicheck (checkfcn, errormsg, errid, varargin)
% performs check on multiple inputs and raises error on failure

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