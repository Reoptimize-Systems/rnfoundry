function settings = mcoreerrormail(settings, errorstruct, i, parameterCell, evalfcn)

    persistent errorid
    
    if isempty(errorid) || ~strcmp(char(errorid), errorstruct.identifier)
        
        errorid = errorstruct.identifier;
        
        % Send email warning about this
        emailstr = sprintf(['This is an automated message from mcore\n\n', ...
            'The objective function has encountered the following error:\n\n%s\n\nWhile evaluating ind %d'], errorstruct.message, i);

        if isa(evalfcn, 'function_handle')
            evalfcn = func2str(evalfcn);
        else
            if ~ischar(evalfcn)
                warning('evalfcn was neither a string nor a function handle, mcore email not being sent.')
                return;
            end
        end

        try
            sendmail('r.crozier@ed.ac.uk', ...
                     sprintf('mcore: error in objective function %s', evalfcn), ...
                     emailstr);
        catch
            warning('mcore objective fcn error email send failure.')
        end
    
    end

end