function struct2tabular(structure, varargin)
% Outputs a latex tabular or threeparttable containing the field names of a
% structure and their contents.
%
% 
    Inputs.fid = 1;
    Inputs.TableType = 'nt';
    Inputs.Caption = '';
    Inputs.CaptionPos = 'top';
    Inputs.Label = '';
    Inputs.NumFormat = '%-6.4g';
    Inputs.MaxMatrixPrintout = 20;
    Inputs.Columns = 2;
    Inputs.ColumnTypes = '';
    Inputs.UseBooktabs = true;
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    if isempty(Inputs.ColumnTypes)
        Inputs.ColumnTypes = repmat('c', 1, Inputs.Columns*2);
    end

    fprintf(Inputs.fid, '\n');
    
    if strcmpi('3pt', Inputs.TableType)
        fprintf(Inputs.fid, '\\begin{threeparttable}\n');
        isthreeparttable = true;
    elseif strcmpi('l3pt', Inputs.TableType)
        fprintf(Inputs.fid, '\\begin{ThreePartTable}\n');
        isthreeparttable = true;
    else
        isthreeparttable = false;
    end
    
    % get the field names of the structure
    fnames = fieldnames(structure);

    tnotecnt = 0;
    
    if strcmpi('l3pt', Inputs.TableType)

        fprintf(Inputs.fid, '\t\\begin{TableNotes}\n');
        
        for i = 1:numel(fnames)

            val = getfield(structure, fnames{i});

            if isnumeric(val) && ~isscalar(val) && (isvector(val) || ismatnotvec(val))
                
                if numel(val) > Inputs.MaxMatrixPrintout
                    fprintf(Inputs.fid, '\t\t\\item[%d] contents: %s \\ldots\n', tnotecnt, sprintf([Inputs.NumFormat, ' '], val));
                else
                    fprintf(Inputs.fid, '\t\t\\item[%d] contents: %s\n', tnotecnt, sprintf([Inputs.NumFormat, ' '], val));
                end
                
                tnotecnt = tnotecnt + 1;
                
            end

        end
        
        fprintf(Inputs.fid, '\t\\end{TableNotes}\n');
    end

    if strcmpi('l3pt', Inputs.TableType)
        fprintf(Inputs.fid, '\t\\begin{longtable}{%s}\n', Inputs.ColumnTypes);
        if ~isempty(Inputs.Caption) && strcmpi(Inputs.CaptionPos, 'top') && ischar(Inputs.Caption)
            fprintf(Inputs.fid, '\\caption{%s}', Inputs.Caption)
        end
        if ~isempty(Inputs.Label) && ischar(Inputs.Label)
            fprintf(1, '\n\\label{%s} \\\\\n', Inputs.Label);
        else
            fprintf(1, '\b\\\\\n');
        end
    else
        if ~isempty(Inputs.Caption) && strcmpi(Inputs.CaptionPos, 'top') && ischar(Inputs.Caption)
            fprintf(Inputs.fid, '\\caption{%s}\n', Inputs.Caption)
        end
        if ~isempty(Inputs.Label) && ischar(Inputs.Label)
            fprintf(1, '\\label{%s}\n', Inputs.Label);
        end
        fprintf(Inputs.fid, '\t\\begin{tabular}{%s}\n', Inputs.ColumnTypes);
    end
    
    if Inputs.UseBooktabs
        fprintf(Inputs.fid, '\t\\toprule\n');
    end
    
    fprintf(Inputs.fid, '\t');
    
    for i = 1:Inputs.Columns
        fprintf(Inputs.fid, 'Field & Contents &');
    end
    
    if strcmpi('l3pt', Inputs.TableType)
        tnotestr = 'tnote';
        if Inputs.UseBooktabs
            fprintf(Inputs.fid, '\b\\\\\n\t\\midrule\n\t\\endhead\n\t\\cmidrule{%d-%d}\n\t\\multicolumn{%d}{r}{\\textit{continued}}\n\t\\endfoot\n\t\\insertTableNotes\\\\\n\t\\endlastfoot\n', Inputs.Columns*2, Inputs.Columns*2, Inputs.Columns*2);
        else
            fprintf(Inputs.fid, '\b\\\\\n\t\\endhead\n\t\\multicolumn{%d}{r}{\\textit{continued}}\n\t\\endfoot\n\t\\insertTableNotes\\\\\n\t\\endlastfoot\n', Inputs.Columns*2);
        end
    else
        tnotestr = 'tnote';
        if Inputs.UseBooktabs
            fprintf(Inputs.fid, '\b \\midrule \\\\\n')
        else
            fprintf(fid, '\b \\\\\n')
        end
    end
    
    matfields = {'dummy'};
    
    tnotecnt = 0;
    
    colcnt = 1;
    
    colstrt = '\t\t';
    
    for i = 1:numel(fnames)

        val = getfield(structure, fnames{i});
        
        fnames{i} = strrep(fnames{i}, '_', '\_');

        if ischar(val)

            fprintf(Inputs.fid, [colstrt, '%s & %s'], fnames{i}, val);

        elseif isnumeric(val)

            if isscalar(val)
                fprintf(Inputs.fid, [colstrt, '%s & ', Inputs.NumFormat], fnames{i}, val);
            else
                if isvector(val)
                    fprintf(Inputs.fid, [colstrt, '%s & vector'], fnames{i});
                    if isthreeparttable
                        matfields{end+1} = fnames{i};
                        tnotecnt = tnotecnt + 1;
                        fprintf(Inputs.fid, '\\%s{%d}', tnotestr, tnotecnt);
                    end
                elseif ismatnotvec(val)
                    fprintf(Inputs.fid, [colstrt, '%s & matrix'], fnames{i});
                    if isthreeparttable
                        matfields{end+1} = fnames{i};
                        tnotecnt = tnotecnt + 1;
                        fprintf(Inputs.fid, '\\%s{%d}', tnotestr, tnotecnt);
                    end
                else
                    fprintf(Inputs.fid, [colstrt, '%s & %s'], fnames{i}, '');
                end
            end

        elseif isstruct(val)
            fprintf(Inputs.fid, [colstrt, '%s & %s'], fnames{i}, 'structure');
        end
    
        if colcnt == Inputs.Columns
            fprintf(Inputs.fid, ' \\\\\n');
            colcnt = 1;
            colstrt = '\t\t';
        else
            
            fprintf(Inputs.fid, ' & ');
            
            if i == numel(fnames)
                rowfinstr = repmat('&', 1, (Inputs.Columns-colcnt)*2);
                fprintf(Inputs.fid, '\b\b%s\\\\\n', rowfinstr);
            else
                
                colcnt = colcnt + 1;
                colstrt = '';
            end
                        
        end
    end
    
    if ~isempty(Inputs.Caption) && strcmpi(Inputs.CaptionPos, 'bottom') && ischar(Inputs.Caption)
        fprintf(Inputs.fid, '\\caption{%s}\n', Inputs.Caption)
    end

    if strcmp('l3pt', Inputs.TableType)
        fprintf(Inputs.fid, '\\end{longtable}\n');
        fprintf(Inputs.fid, '\\end{ThreePartTable}\n');
    else
        fprintf(Inputs.fid, '\t\\end{tabular}\n');
    end

    if strcmp('3pt', Inputs.TableType)
        
        fprintf(Inputs.fid, '\t\\begin{tablenotes}\n');

        for i = 1:tnotecnt

            val = getfield(structure, matfields{i+1});

            if numel(val) > Inputs.MaxMatrixPrintout
                fprintf(Inputs.fid, '\t\t\\item[%d] contents: %s \\ldots\n', tnotecnt, sprintf([Inputs.NumFormat, ' '], val));
            else
                fprintf(Inputs.fid, '\t\t\\item[%d] contents: %s\n', tnotecnt, sprintf([Inputs.NumFormat, ' '], val));
            end

            fprintf(Inputs.fid, '\t\t\\item[%d] %s\n', i, notetext);

        end
       
       fprintf(Inputs.fid, '\t\\end{tablenotes)\n\\end{threeparttable}\n');
    end
            
        
% \begin{threeparttable}
% \caption{} \label{putalablehere}
% \begin{tabular}
% Table contents here
% you can add tnote{putasignhere} for some thing
% \end{tabular}
% \begin{tablenotes}
% \item[]
% \item[putasign_here]
% \end{tablenotes}
% \end{threeparttable}

end