function outCell = mcoresimfun_ACTM(design, simoptions, options)
% mcoresimfun_ACTM: calls simfun_ACTM and returns the results in a cell
% array suitable to be passed to startmulticoremaster2.m

    [design, simoptions] = simfun_ACTM(design, simoptions);

    outCell = {design, simoptions, options};
    
end