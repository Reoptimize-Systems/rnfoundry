function dispstruct(s, colwid, tabwid)
% displays the contents of a structure in a nicely formatted fashion
%
% Syntax
%
% dispstruct(s, colwid, tabwid)
%
% Description
%
% dispstruct displays the contents of a strucuture at the command line in a
% readable format. Only scalar structures can be displayed.
%
% First all fields containing nonscalar numeric values or strings are
% displayed as short one-line descriptions of their type and possibly a few
% values of their contents, depending on their type. Next the scalar
% members and strings are displayed as a multi-column table.
%

    if nargin < 2
        colwid = 10;
    end
    
    if nargin < 3
        tabwid = 80;
    end
    
    if ~isscalar(s)
        error('only scalar structure arrays are supported')
    end
    
    fnames = fieldnames(s);
    
    % construct the table contents, or display info on fields
    hasnonscalar = false;
    
    fprintf(1, 'Nonscalars and non character arrays:\n');
    
    tablefnames = {};
    tablecontents = {};
    
    for i = 1:numel(fnames)
       
        if (isnumeric(s.(fnames{i})) && isscalar(s.(fnames{i}))) || ischar(s.(fnames{i}))
           
            tablefnames = [tablefnames; fnames(i)];
            
            tablecontents = [tablecontents; {s.(fnames{i})}];
            
        else
            
            hasnonscalar = true;

            % display the info
            fprintf(1, 'Field: %s\t%s\t%s\t', fnames{i}, mat2str(size(s.(fnames{i}))), class(s.(fnames{i})))
            
            if isnumeric(s.(fnames{i}))    
                
                dispvals = s.(fnames{i})(1:min(4,numel(s.(fnames{i}))));
                
                fprintf(1, '[ ');
                
                for j = 1:numel(dispvals)
                    fprintf(1, '%g, ', dispvals(j));
                end
                
                fprintf(1, '... ]\n');
            else
                fprintf(1, '\n');
            end
            
        end
        
    end
    
    if ~hasnonscalar
       fprintf(1, 'No nonscalars or non character arrays\n') 
    end
    
    if numel(tablefnames) > 0
        
        % choose a suitible number of table columns
        tabcols = 2*max(1, round((tabwid/2) / colwid));

        if numel(tablefnames) < (tabcols/2)
            tabcols = 2 * numel(tablefnames);
        end
    
        % construct the displaytable arguments appropriately
        
        fprintf(1, 'Scalar fields and character arrays:\n')
    
        tabrows = ceil(numel(tablefnames) / (tabcols/2));
        
        M = mod(numel(tablefnames), (tabcols/2));
        
        if M ~= 0
            M = (tabcols/2) - M;
        end
        
        tablefnames = [tablefnames; cell(M, 1)];
        tablecontents = [tablecontents; cell(M, 1)];
        
        % prepare the table data cell array
        data = cell(tabrows, tabcols);
        
        % construct the column headings
        colheadings = repmat({'Field', 'Contents'}, 1, tabcols/2);

        % put the field names into the table data
        data(:,1:2:end-1) = reshape(tablefnames,  tabcols/2, tabrows)';
        
        % put the values into the table data
        data(:,2:2:end) = reshape(tablecontents, tabcols/2, tabrows)';
        
        % display the table
        displaytable(data, colheadings, colwid)
        
    else
        
        fprintf(1, 'No scalars or character arrays\n') 

    end

end