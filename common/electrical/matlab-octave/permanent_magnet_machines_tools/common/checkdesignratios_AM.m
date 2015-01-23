function design = checkdesignratios_AM (design, ratiospec, dumpstruct)
% checks design ratios are within specified limits

    if nargin < 3
        dumpstruct = false;
    end
    
    emax = @(fn) fprintf (1, 'Ratio field %s exceeded the specified maximum possible allowed\n', fn);
    emin = @(fn) fprintf (1, 'Ratio field %s is below the specified minimum possible allowed\n', fn);
    
    results = checkstructratios (design, ratiospec, dumpstruct);
    
    for ind = 1:size (ratiospec, 1)
        
       if results(ind) == 1
           emin (ratiospec{ind,1});
       elseif results(ind) == 2
           emax (ratiospec{ind,1});
       elseif results(ind) == 3
           emin (ratiospec{ind,1});
           emax (ratiospec{ind,1});
       end
       
    end
    
    if any (results ~= 0)
        
        if dumpstruct
            dispstruct (design);
        end
        
        error ('RENEWNET:pmmachines:designspec', ...
               'Impossible design spec');
    end
    
end


function result = checkstructratios (S, ratiospec)

    result = zeros (size (ratiospec, 1), 1);
    
    for ind = 1:size (ratiospec, 1)
        
        if isfield (S, ratiospec{ind,1})
            
            if  S.(ratiospec{ind,1}) < ratiospec{ind,2}
                
                result(ind) = result(ind) + 1;
                
                fprintf (1, 'Struct Ratio %s was %f, which was less than the specified min possible allowed (%f)\n', ...
                            ratiospec{ind,1}, ratiospec{ind,2}, S.(ratiospec{ind,1}));
                
            end
            
            if  S.(ratiospec{ind,1}) > ratiospec{ind,3}
                
                result(ind) = result(ind) + 2;
                
                fprintf (1, 'Struct Ratio %s was %f, which was less than the specified min possible allowed (%f)\n', ...
                            ratiospec{ind,1}, ratiospec{ind,3}, S.(ratiospec{ind,1}));
                
            end
            
        else
            result(ind) = 4;
        end  
        
    end

end