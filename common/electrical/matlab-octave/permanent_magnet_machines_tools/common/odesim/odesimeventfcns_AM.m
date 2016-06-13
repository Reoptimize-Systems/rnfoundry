function [value,isterminal,direction] = odesimeventfcns_AM (t, y, design, simoptions)
% evaluate a collection of outputfcns in an ode simulation
    
    compnames = fieldnames (simoptions.ODESim.SolutionComponents);

    value = [];
    isterminal = [];
    direction = [];
    
    for ind = 1:numel (compnames)

        if isfield (simoptions.ODESim.SolutionComponents.(compnames{ind}), 'EventFcn') ...
                && ~isempty (simoptions.ODESim.SolutionComponents.(compnames{ind}).EventFcn)
            
                [tmpvalue,tmpisterminal,tmpdirection] = ...
                    feval ( simoptions.ODESim.SolutionComponents.(compnames{ind}).OutputFcn, ...
                            t, ...
                            y(simoptions.ODESim.SolutionComponents.(compnames{ind}).SolutionIndices,end), ...
                            design, ...
                            simoptions );
                             
                value = [ value, tmpvalue];
                isterminal = [ isterminal, tmpisterminal ];
                direction = [ direction, tmpdirection];

        end

    end

end