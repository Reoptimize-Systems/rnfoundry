function outCell = mcoresimfun_AM(design, simoptions, varargin)
% mcoresimfun_AM: calls simoptions.ODESim.PreProcFcn and returns the results in a cell
% array suitable to be returned to startmulticoremaster2.m


    simoptions = rmiffield(simoptions, 'filenamebase');
    
    [design, simoptions] = feval(simoptions.ODESim.PreProcFcn, design, simoptions);
    
    simoptions.ODESim.PreProcFcn = [];
    
    outCell = [{design, simoptions}, varargin];
    
end