function outCell = mcoresimfun_AM(design, simoptions, varargin)
% mcoresimfun_AM: calls simoptions.simfun and returns the results in a cell
% array suitable to be returned to startmulticoremaster2.m


    simoptions = rmiffield(simoptions, 'filenamebase');
    
    [design, simoptions] = feval(simoptions.simfun, design, simoptions);
    
    simoptions.simfun = [];
    
    outCell = [{design, simoptions}, varargin];
    
end