function s = cellstr2str(C)
% converts a cellarray of strings to a single string with newline
% characters between each string

    s = sprintf('%s\n', C{:});

end