function outCell = mcoresimfun_PMSM(design, simoptions, options)
% mcoresimfun_ACPMSM: calls simfun_ACPMSM and returns the results in a cell
% array suitable to be passed to startmulticoremaster2.m

    [design, simoptions] = simfunnocurrent_PMSM(design, simoptions);

    outCell = {design, simoptions, options};
    
end