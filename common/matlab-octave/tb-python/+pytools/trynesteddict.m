function cellstrs = trynesteddict (dictname, nestednames, assignment)

nname = pytools.nesteddictname (nestednames);

cellstrs = { ...
        'try:';
sprintf('    %s%s = %s', dictname, nname, assignment);
        'except KeyError:';
sprintf('    %s%s = {}', dictname, py.nesteddictname (nestednames(1:end-1)));
sprintf('    %s%s = %s', dictname, nname, assignment);
};

end
