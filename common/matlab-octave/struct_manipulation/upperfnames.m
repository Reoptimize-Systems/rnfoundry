function s = upperfnames(s)
% upperfnames: converts all the field names in a structure to upper case

    % get the structure field names
    fn = fieldnames(s);

    for ii = 1:length(fn)

        fname = fn{ii};
        
        % check if the field is already upper case
        if ~strcmp(fname, upper(fname))
            
            if isfield(s, upper(fname))
                error('Both upper case and lower case or mixed case versions of the same field name are already present in structure.')
            else
                % if not create an upper case copy of the field
                s = RenameField(s, fname, upper(fname));
            end
        
        end

    end

end