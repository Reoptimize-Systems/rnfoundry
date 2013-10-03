function outCell = mcoresimfun_ACTIAM(design, simoptions, options)
% mcoresimfun_ACTM: calls simfun_ACTM and returns the results in a cell
% array suitable to be passed to startmulticoremaster2.m

    [design, simoptions] = simfun_ACTIAM(design, simoptions);

    outCell = {design, simoptions, options};
    
end