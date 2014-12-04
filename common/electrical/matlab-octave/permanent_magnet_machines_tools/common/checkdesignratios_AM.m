function design = checkdesignratios_AM (design, ratiospec)
% checks design ratios are within specified limits
    
    emax = @(fn) fprintf (1, 'Ratio field %s exceeded the specified maximum possible allowed\n', fn);
    emin = @(fn) fprintf (1, 'Ratio field %s is below the specified minimum possible allowed\n', fn);
    
    results = checkstructratios (design, ratiospec);
    
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
                
            end
            
            if  S.(ratiospec{ind,1}) > ratiospec{ind,3}
                
                result(ind) = result(ind) + 2;
                
            end
            
        else
            result(ind) = 4;
        end  
        
    end

end