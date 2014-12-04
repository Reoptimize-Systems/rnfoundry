function pystring = nesteddictname (nestednames)

    if ischar (nestednames)
        % convert to a cell array of strings
        nestednames = {nestednames};
    end
    
    pystring = '';
    for ind = 1:numel(nestednames)
       
        pystring = [pystring, sprintf('["%s"]', nestednames{ind})];
        
    end
    
end